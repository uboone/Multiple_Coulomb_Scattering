#pragma once

#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TMatrixD.h"
#include "TVector3.h"
#include "TVectorT.h"
#include "TArrayD.h"
#include "TMatrixDEigen.h"

struct ComparePCAProjection;
struct Comparator;
struct Point;
struct Track;

namespace MCSHelper {

  std::vector<double> diff (std::vector<double> u, std::vector<double> v) {
    for (int i=0,n=u.size();i<n;i++) { u[i] -= v[i]; }
    return u;
  }

  double norm (std::vector<double> v) {
    double norm = 0;
    for (int i=0,n=v.size();i<n;i++) { norm += std::pow(v[i],2); }
    return std::sqrt(norm);
  }

  double dot (std::vector<double> u, std::vector<double> v) {
    double val = 0;
    for (int i=0,n=u.size();i<n;i++) { val += u[i]*v[i]; }
    return val;
  }

  std::vector<double> cross (std::vector<double> u, std::vector<double> v) {
    return { u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0] };
  }

  std::vector<double> scale (std::vector<double> v, double a) {
    for (int i=0,n=v.size();i<n;i++) { v[i] *= a; }
    return v;
  }

  //compares two point-distance tuples based on their pre-computed distance from a reference point
  bool compare_edges (std::tuple<double,Point*> t1, std::tuple<double,Point*> t2) { return std::get<0>(t1) < std::get<0>(t2); }

}

//wrapper structure that lets you sort points with respect to an axis
struct ComparePCAProjection {  
  std::vector<double> axis;
  ComparePCAProjection() {
    (*this).axis = {0,0,0};
  }

  //Helper function for sorting points
  bool compare(std::vector<double> a, std::vector<double> b){ return MCSHelper::dot(a,(*this).axis)<MCSHelper::dot(b,(*this).axis); }

  //Sort track points in terms of PCA axis
  void sort_points(std::vector<std::vector<double>> &points, std::vector<double> principleAxis){
    (*this).axis = principleAxis;
    //std::sort(points.begin(), points.end(), (*this).compare);
    std::sort(points.begin(), points.end(), [this](std::vector<double> a, std::vector<double> b) { return MCSHelper::dot(a,(*this).axis)<MCSHelper::dot(b,(*this).axis); });
  }
};

//for use in forming "shortest" path along trajectory points
struct Point {
  int id;
  bool exhausted = false;
  double score;
  std::vector<double> pos, edges_dist;
  Point* prior = NULL;
  std::vector<Point*> edges;

  Point (int id, std::vector<double> pos) {
    (*this).id = id;
    (*this).pos = pos;
    (*this).prior = NULL;
    (*this).score = 1e9;
  }

  //assigns a score to a distance (in cm).  The linear term encourages a short path while the quadratic term demands small distances between segments.
  //These exponents could be changed, particularly if 0.6 cm spacing between points is not used.
  double get_dist_score (Point p) {
    double dist = MCSHelper::norm(MCSHelper::diff((*this).pos,p.pos));
    return dist + std::pow(dist,2);
  }

  //attempts to update the path to point p if the route through the current point is better than the previous route
  void update_path (int edge_index) {
    Point* edge = (*this).edges[edge_index];
    double new_score = (*this).score + (*this).edges_dist[edge_index];
    if (new_score < (*edge).score) {
      (*edge).score = new_score;
      (*edge).prior = this;
      (*edge).exhausted = false;
    }
  }

  //attempts to update the path to all edges
  void update_paths () {
    for (int i=0,n=(*this).edges.size();i<n;i++) { update_path(i); }
    (*this).exhausted = true;
  }

  void sort_edges () {
    std::vector<std::tuple<double,Point*>> edges_tuples;
    for (int i=0,n=(*this).edges.size();i<n;i++) { edges_tuples.push_back( std::make_tuple((*this).edges_dist[i],(*this).edges[i]) ); }
    std::sort(edges_tuples.begin(), edges_tuples.end(), MCSHelper::compare_edges);
    for (int i=0,n=edges_tuples.size();i<n;i++) {
      (*this).edges_dist[i] = std::get<0>(edges_tuples[i]);
      (*this).edges[i]      = std::get<1>(edges_tuples[i]);
    }
  }

  //only add edges along the muon direction (reco end-start) and only keep 4 closest edges
  void add_edge (Point* p, std::vector<double> dir) {
    int nedges_max = 20;
    bool in_dir = MCSHelper::dot(MCSHelper::diff((*p).pos,(*this).pos), dir) > 0;
    double dist = get_dist_score(*p);
    int nedges = (*this).edges_dist.size();
    if (in_dir && (nedges<nedges_max || dist<(*this).edges_dist.back())) {
      (*this).edges.push_back(p);
      (*this).edges_dist.push_back(dist);
      sort_edges();
      if (nedges>=nedges_max) {
        (*this).edges.pop_back();
        (*this).edges_dist.pop_back();
      }
    }
  }

};

//compares two points based on their score (for use in forming "shortest" path)
struct Comparator {
  bool operator()(const Point* p1, const Point* p2) const { return !(*p1).exhausted && ((*p2).exhausted || ((*p1).score < (*p2).score)); }
};

struct Track {
  std::vector<std::vector<double>> points;
  std::vector<double> weights;
  double total_weight;
  int N; // number of vertices 

  Track () {
    (*this).N            = 0;
    (*this).total_weight = 0;
  }

  void add_point (std::vector<double> point, double weight=1.) {
    (*this).points.push_back(point);
    (*this).weights.push_back(weight);
    (*this).total_weight += weight;
    (*this).N++;
  }

  //removes segment from a track.  Assumes segment and track are sorted the same
  void remove_seg (std::vector<double> first_pos, int size) {
    //find first pos
    int first_index = -1;
    for (int i=0;i<(*this).N;i++) {
      std::vector<double> point_pos = (*this).points[i];
      if (point_pos[0]!=first_pos[0] || point_pos[1]!=first_pos[1] || point_pos[2]!=first_pos[2]) { continue; }
      first_index = i;
      break;
    }

    remove_seg_at(first_index,size);
  }

  //remove a segment of specified size at a specified index
  void remove_seg_at (int first_index, int size) {
    double removed_weight = 0;
    for (int i=0;i<size;i++) { removed_weight += (*this).weights[first_index+i]; } 
    (*this).points.erase( (*this).points.begin()+first_index,  (*this).points.begin() +first_index+size);
    (*this).weights.erase((*this).weights.begin()+first_index, (*this).weights.begin()+first_index+size);
    (*this).total_weight -= removed_weight;
    (*this).N            -= size;

  }

  void clear () {
    (*this).points.clear();
    (*this).weights.clear();
    (*this).total_weight = 0;
    (*this).N = 0;
  }
};

//==============================================================================================================================================================

class MCS{
  public:

    MCS();
    void cleanUp();
    void run(std::vector<double> vtx_muon_start_reco, std::vector<double> vtx_muon_end_reco, std::vector<std::vector<double>> trajectory_points_initial);

  protected:

    //helper functions
    void configure();
    void writeOutput();
    double get_angle(std::vector<double> v1, std::vector<double> v2);
    void setSegAngles(std::vector<std::vector<double>> priorSegFitVector, std::vector<double> segFitVector, std::vector<double> &angle_vec, std::vector<double> &angleProjX_vec, std::vector<double> &angleProjY_vec);
    std::vector<std::vector<double>> fitPCA(Track track,  std::vector<double> &com, std::vector<double> &evals);
    Track get_seg (Track track, std::vector<double> axis, double seg_length);
    bool fitSegPCA(Track &track, std::vector<double> aAxis, Track &seg, std::vector<std::vector<double>> &segFitVector, std::vector<double> &segCOM, double seg_length, double &currentFirstPointProj, double &currentLastPointProj);
    std::tuple<bool,std::vector<std::vector<double>>> trim_trajectory(double npoints_traj, std::vector<std::vector<double>> trajectory_points_initial, std::vector<double> vtx_muon_start_reco, std::vector<double> vtx_muon_end_reco);
    std::tuple< std::vector<Track>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> > form_segs (std::vector<std::vector<double>> vec_points, std::vector<double> muon_start, std::vector<double> muon_end, double seg_length);

    double beta(double gamma);
    double gamma(double KE, double mass);
    double sigmoid(double x);
    void increment_energy (double &e);
    void decrement_dist (double &x, double xmax);
    void setUKEfromRR();
    void setUKEfromEX (TGraph *uKEfromRR, TGraph *uRRfromKE);
    double sigmaH (double T);
    double sigmoid(double x, std::vector<double> par);
    double quartic_decay(double x, std::vector<double> par);
    std::vector<double> pred_theta_xz_pars(double T);
    std::vector<double> pred_theta_yz_pars(double T, int vx_index);
    double double_gaussian(double angle, std::vector<double> pars);
    double lnlikelihood_theta_xz(double angle, double T);
    double lnlikelihood_theta_yz(double angle, double T, double vx);
    double lnlikelihood_track(double* KE, double* par);
    std::vector<double> estimate_energy(std::vector<double> segs_distance, std::vector<double> segs_angle_x, std::vector<double> segs_angle_y,  std::vector<double> vx_comps);

    //constants
    const double seg_length = 14; //cm
    const double I = 188.0*pow(10,-6); 			//eV
    const double Mmu = 105.658;				// MeV for muon
    const double rho = 1.396;				// LAr density [g/cm3]

    //variables set in fhicl file
    bool f_wirecellPF;
    std::string fPFInputTag;
    std::string fportedWCSpacePointsTrecchargeblobLabel;
    double res_sigma1_xz, res_sigma2_xz;
    std::vector<double> res_sigma1_yz, res_sigma2_yz, par_sigma1_xz, par_sigma2_xz, par_ratio_xz;
    std::vector<std::vector<double>> par_sigma1_yz, par_sigma2_yz, par_ratio_yz;

    //internal variables
    TGraph* uKEfromRR;
    TGraph* uRRfromKE;
    double mu_tracklen;
    double emu_tracklen;
    double emu_MCS;
    double ambiguity_MCS;
};


