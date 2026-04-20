#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include <memory>
#include <string>
#include <dirent.h>
//#include <iostream>
#include <numeric>
#include "TGraph2D.h"
#include "TGraph.h"
#include "TF1.h"

class MCS;

//Constructor
MCS::MCS()
{
  //Create graphs mapping between residual range and energy
  (*this).setUKEfromRR();

  (*this).configure();
}

//initialization using fhicl parameters
void MCS::configure()
{
  //Double Gaussian tune parameters
  res_sigma1_xz = 0.005776;
  res_sigma2_xz = 0.01821;
  par_sigma1_xz = { -0.449144931, 0.793132642, -1.291292240, 0.536765147, -0.084910516, 0.146304242 };
  par_sigma2_xz = {  0.562850793, 0.118940108,  0.000625509, 0,           -0.000000100, 1.217639251 };
  par_ratio_xz  = {  0.839684805, 0.839684805,  0, 0 };
  res_sigma1_yz = { 0.0449, 0.0206,  0.01403, 0.0131,  0.0114  };
  res_sigma2_yz = { 0.1506, 0.06154, 0.03965, 0.04179, 0.07347 };
  par_sigma1_yz = {{ -0.09,        -0.084325217,  0.487240052,  0.395496655, -0.187184874, 0.166128734 },
                   {  0.0,         -0.575280374,  0.070151974,  0.187260875, -0.099717108, 0.160128002 },
                   { -0.153367057, -0.583042532,  0.983374136, -0.712652874,  0.134743902, 0.465439107 },
                   { -0.268993212,  0.103899779, -0.588953942,  0.282356661, -0.067930741, 0.176282668 },
                   { -0.2,         -0.724028910,  0.660065851, -0.327141529,  0.038426745, 0.094357571 }};
  par_sigma2_yz = {{  0.193250363,  3.046073125,  15.0,        -7.519451962, -1.0,         0.5         },
                   {  2.226214634, -5.407245785,  7.227026554, -2.411882646, -0.000001517, 0.224728907 },
                   {  0.25,         1.670060693, -2.086639043,  1.122640876,  0.000000440, 0.268340031 },
                   {  0.430965414, -0.079204912,  0.220949488, -0.000010000, -0.000010000, 0.3         },
                   {  1.635001641, -2.414781711,  1.669221724, -0.320484259,  0.0,         0.157018001 }};
  par_ratio_yz  = {{  0.519230936,  1.0,          1.663982350,  0.084218766 },
                   {  0.7,          0.788705752,  5.0,          2.045055158 },
                   {  0.822433461,  0.701381849,  5.0,          0.0         },
                   {  0.897408009,  0.7,          2.182533745,  0.297834161 },
                   {  0.9,          0.9,          0.0,          0.0         }};
}


//clean up allocated memory
void MCS::cleanUp()
{
  uKEfromRR = nullptr;
  uRRfromKE = nullptr;
}

//Prints out MCS energy estimate and ambiguity score (1=most ambiguous)
void MCS::writeOutput(){
  std::cout << "emu_MCS, ambiguity_MCS = " << emu_MCS << ",  " << ambiguity_MCS << std::endl;
}

void MCS::run(std::vector<double> vtx_muon_start_reco, std::vector<double> vtx_muon_end_reco, std::vector<std::vector<double>> trajectory_points_initial){

  std::cout << "Begin MCS::run" << std::endl;

  //prepare outputs
  mu_tracklen   = -1;
  emu_tracklen  = -1;
  emu_MCS       = -1;
  ambiguity_MCS = -1;

  //Trim to only include muon path (no delta rays or crossing tracks)
  int npoints_trajectory = trajectory_points_initial.size();
  std::tuple<bool,std::vector<std::vector<double>>> trajectory_tuple = MCS::trim_trajectory(npoints_trajectory, trajectory_points_initial, vtx_muon_start_reco, vtx_muon_end_reco);
  bool bad_path                                            = std::get<0>(trajectory_tuple);
  std::vector<std::vector<double>> trajectory_points_final = std::get<1>(trajectory_tuple);
  int npoints_trajectory_final = trajectory_points_final.size();

  //compute residual range and emu_rr
  double rr_path = 0;
  for (int i=1;i<npoints_trajectory_final;i++) { rr_path += MCSHelper::norm(MCSHelper::diff(trajectory_points_final[i],trajectory_points_final[i-1])); }
  double KE_rr_path = uKEfromRR->Eval(rr_path);
  mu_tracklen = rr_path;
  emu_tracklen = (KE_rr_path+Mmu)/1000.;   //convert from KE to E and from MeV to GeV

  //skip events where a path cannot be traversed from muon start to end
  //skip events with very short tracks (remember the track is in reverse order)
  if (bad_path || npoints_trajectory_final<20 || MCSHelper::norm(MCSHelper::diff(trajectory_points_final.back(),vtx_muon_end_reco)) < 2*seg_length) {
    std::cout << "Bad/short path" << std::endl;
    writeOutput();
    return;
  }

  //segment path
  std::cout << "Segmenting muon path" << std::endl;
  std::tuple< std::vector<Track>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> > segs_tuple = MCS::form_segs(trajectory_points_final, vtx_muon_start_reco, vtx_muon_end_reco, seg_length);             
  std::vector<Track> segs                     = std::get<0>(segs_tuple);
  std::vector<std::vector<double>> axes       = std::get<1>(segs_tuple);
  std::vector<std::vector<double>> segs_aAxes = std::get<2>(segs_tuple);
  std::vector<std::vector<double>> segs_COM   = std::get<3>(segs_tuple);
  std::vector<double> segs_distance           = std::get<4>(segs_tuple);  //distance from muon start to seg midpoint
  std::vector<double> segs_angle              = std::get<5>(segs_tuple);
  std::vector<double> segs_angle_projB        = std::get<6>(segs_tuple);
  std::vector<double> segs_angle_projC        = std::get<7>(segs_tuple);
  //std::vector<double> segs_displacement       = std::get<8>(segs_tuple);

  //skip events without any angle measurements
  if (segs.size()<2) {
    std::cout << "Short path" << std::endl;
    writeOutput();
    return;
  }

  //Get the x-component of each segment fit vector
  std::vector<double> vx_components;
  for (const auto& seg_dir : segs_aAxes) {
    if (!seg_dir.empty()) { vx_components.push_back(seg_dir[0]); }
  }

  //estimate energy
  std::cout << "Estimating energy" << std::endl;
  std::vector<double> kemu_MCS_tuple = estimate_energy(segs_distance, segs_angle_projB, segs_angle_projC, vx_components);
  emu_MCS       = (kemu_MCS_tuple[0]+Mmu)/1000.;   //convert from KE to E and from MeV to GeV
  ambiguity_MCS =  kemu_MCS_tuple[1];

  //write outputs
  writeOutput();
  std::cout << "Done MCS" << std::endl;
}


//Find the angle separation between two vectors in 3d
double MCS::get_angle(std::vector<double> v1, std::vector<double> v2) {
  double l1 = MCSHelper::norm(v1);
  double l2 = MCSHelper::norm(v2);
  return (acos( MCSHelper::dot(v1,v2)/(l1*l2)) );
}

//Helper function that computes the angle between pca fit segments and stores the 3D and 2D angles
void MCS::setSegAngles(std::vector<std::vector<double>> priorSegFitVector, std::vector<double> segFitVector,
                  std::vector<double> &angle_vec, std::vector<double> &angleProjX_vec, std::vector<double> &angleProjY_vec){

  std::vector<double> aAxis_prior = priorSegFitVector[0];
  std::vector<double> bAxis_prior = priorSegFitVector[1];
  std::vector<double> cAxis_prior = priorSegFitVector[2];


  std::vector<double> vecy_plane  = MCSHelper::cross(aAxis_prior,{1,0,0});
  std::vector<double> vecx_plane  = MCSHelper::cross(aAxis_prior,vecy_plane);
  vecx_plane = MCSHelper::scale(vecx_plane,1./MCSHelper::norm(vecx_plane));
  vecy_plane = MCSHelper::scale(vecy_plane,1./MCSHelper::norm(vecy_plane));
  std::vector<double> projX = MCSHelper::diff(segFitVector, MCSHelper::scale(vecy_plane,MCSHelper::dot(segFitVector,vecy_plane))) ;
  std::vector<double> projY = MCSHelper::diff(segFitVector, MCSHelper::scale(vecx_plane,MCSHelper::dot(segFitVector,vecx_plane)));
  
  projX = MCSHelper::scale(projX,1./MCSHelper::norm(projX));
 
  int dirX = 1 - 2*(MCSHelper::dot(segFitVector,vecx_plane)<0);
  int dirY = 1 - 2*(MCSHelper::dot(segFitVector,vecy_plane)<0);

  angle_vec.push_back(MCS::get_angle(segFitVector,aAxis_prior));
  angleProjX_vec.push_back( MCS::get_angle(projX, aAxis_prior)*dirX );
  angleProjY_vec.push_back( MCS::get_angle(projY, aAxis_prior)*dirY );
}

//Track t is set of 3D points, m is to be computed as mean value of points, the remainder is fit using PCA and stored in a,b,c
std::vector<std::vector<double>> MCS::fitPCA(Track track,  std::vector<double> &com, std::vector<double> &evals){
    
  std::vector<std::vector<double>> meanPoints(track.N);
  TArrayD covData(9);
  TMatrixD* covMatrix = new TMatrixD(3,3);
  TVectorD eigenValues;
  evals.resize(3);
	
  com = {0,0,0};
  for(int i=0; i<track.N; i++){  //compute center of mass
    com[0] += track.points[i][0]*track.weights[i];
    com[1] += track.points[i][1]*track.weights[i];
    com[2] += track.points[i][2]*track.weights[i];
  }
  com = MCSHelper::scale(com,1./track.total_weight);
  
  for(int i=0; i<track.N; i++){	//create mean value matrix
    meanPoints[i].resize(3);
    meanPoints[i][0] = (track.points[i][0]-com[0]) * track.weights[i];
    meanPoints[i][1] = (track.points[i][1]-com[1]) * track.weights[i];
    meanPoints[i][2] = (track.points[i][2]-com[2]) * track.weights[i];		
  }
	
  //compute and put covariance data into a 1D array
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<track.N; k++) {covData[3*i+j] += meanPoints[k][i] * meanPoints[k][j];}
      covData[3*i+j] = covData[3*i+j] / track.N;
    }
  }
  covMatrix->SetMatrixArray(covData.GetArray());
  
  //find principle axis
  std::vector<std::vector<double>> axes = {{0,0,0},{0,0,0},{0,0,0}};
  const TMatrixD eigenVectors = covMatrix->EigenVectors(eigenValues);

  for (int i=0;i<3;i++){
    axes[0][i] = eigenVectors[i][0];
    axes[1][i] = eigenVectors[i][1];
    axes[2][i] = eigenVectors[i][2];
    evals[i] = eigenValues[i];
  }
  evals = MCSHelper::scale(evals, 1./MCSHelper::norm(evals));
  delete covMatrix;
  return axes;
}

//helper that takes a track (set of points) and selects a subset of those points to form a track segment, based on the input axis and length.
Track MCS::get_seg (Track track, std::vector<double> axis, double seg_length) {
  Track track_seg = Track();
	
  double first_proj = MCSHelper::dot(track.points.front(),axis);
  for(int i=0; i<track.N; i++){
    double proj = MCSHelper::dot(track.points[i],axis);
    if (proj < (first_proj+seg_length)){
      track_seg.add_point(track.points[i]);
    }
  }

  return track_seg;
}

//Helper Function that fits a segment with all the points in the next seg_length=14 cm starting after the prior lastPointProj
//points is a vector of all the 3D track points
//aAxis is the primary pca axis vector
//seg is a track to be filled with the current selection of points
//segFitVector will store the value of the pca fit for the seg
//segCOM will store the center of mass for the seg
//PointProj variables store info on where to start and end the segments
bool MCS::fitSegPCA(Track &track, std::vector<double> aAxis, Track &seg, std::vector<std::vector<double>> &segFitVector, std::vector<double> &segCOM,
               double seg_length, double &currentFirstPointProj, double &currentLastPointProj){
  //setup
  seg.points.clear();
  seg.N = 0;
	
  //track_buffer stores only a short amount of the track to speed up computation when sorting
  Track track_buffer = MCS::get_seg(track, aAxis, seg_length*2);

  //update currentFirstPointProj
  double lastPointProj = MCSHelper::dot(track_buffer.points.back(),aAxis);
  currentFirstPointProj = lastPointProj;
  for(int i=0; i<int(track_buffer.points.size()); i++){
    double proj = MCSHelper::dot(track_buffer.points[i],aAxis);
    if((proj>currentLastPointProj) && (proj<currentFirstPointProj)){ currentFirstPointProj = proj; }
  }

  //fill seg with points
  seg = MCS::get_seg(track_buffer, aAxis, seg_length);
  if (seg.N<=1) { return false; }

  //run PCA on the segment and orient along aAxis
  std::vector<double> eigenvalues;
  segFitVector = MCS::fitPCA(seg, segCOM, eigenvalues);
  if(MCSHelper::dot(segFitVector[0],aAxis)<0) {
    segFitVector[0][0] *= -1;
    segFitVector[0][1] *= -1;
    segFitVector[0][2] *= -1;	
  }

  //Sort points along new fit axis
  ComparePCAProjection comparator = ComparePCAProjection();
  comparator.sort_points(track_buffer.points, segFitVector[0]);

  //re-determine what points belong in the segment based on new u-vector
  seg = MCS::get_seg(track_buffer, segFitVector[0], seg_length);
  if (seg.N<=1) { return false; }

  //refit and orient along aAxis
  segFitVector = MCS::fitPCA(seg, segCOM, eigenvalues);
  if(MCSHelper::dot(segFitVector[0],aAxis)<0) {
    segFitVector[0][0] *= -1;
    segFitVector[0][1] *= -1;
    segFitVector[0][2] *= -1;	
  }

  //update current-last-point-proj
  comparator.sort_points(seg.points, segFitVector[0]);
  currentLastPointProj = MCSHelper::dot(seg.points.back(),aAxis);

  //remove points from original track
  comparator.sort_points(seg.points, aAxis);
  track.remove_seg(seg.points.front(),seg.N);
  return true;
}


//helper function that trims an input list of trajectory points to only include those along the shortest path from muon start to end vertices
std::tuple<bool,std::vector<std::vector<double>>> MCS::trim_trajectory(double npoints_traj, std::vector<std::vector<double>> trajectory_points_initial, std::vector<double> vtx_muon_start_reco, std::vector<double> vtx_muon_end_reco) {

  //trim trajectory points to get shortest path along muon track
  //find starting and ending trajectory points based on proximinty to reco muon start and end 
  std::vector<double> muon_dir_reco = MCSHelper::diff(vtx_muon_end_reco,vtx_muon_start_reco);
  int startpoint_index = -1;
  int endpoint_index = -1;
  double startpoint_dist = 1e9;
  double endpoint_dist = 1e9;
  for (int i=0;i<npoints_traj;i++) {
    std::vector<double> pos = trajectory_points_initial[i];
    double dist_start = std::sqrt( MCSHelper::norm(MCSHelper::diff(pos,vtx_muon_start_reco)) );
    double dist_end   = std::sqrt( MCSHelper::norm(MCSHelper::diff(pos,vtx_muon_end_reco)) );
    if (dist_start < startpoint_dist) {
      startpoint_index = i;
      startpoint_dist = dist_start;
    }
    if (dist_end < endpoint_dist) {
      endpoint_index = i;
      endpoint_dist = dist_end;
    }
  }
  //create a vector of all trajectory points, with the 0th Point as the first in the vector. Give the starting point a graph score of 0.
  std::vector<Point*> traj_points;
  Point* new_point0 = new Point(startpoint_index+1, { trajectory_points_initial[startpoint_index][0], trajectory_points_initial[startpoint_index][1], trajectory_points_initial[startpoint_index][2] });
  traj_points.push_back(new_point0);
  (*traj_points[0]).score = 0;
  for (int i=0;i<npoints_traj;i++) {
    Point* new_point = new Point(i+1, { trajectory_points_initial[i][0], trajectory_points_initial[i][1], trajectory_points_initial[i][2] });
    if (i != startpoint_index) { traj_points.push_back(new_point); }
  }
  //form graph of all Points by taking each Point i and calling add_edge on each other point j
  for (int i=0;i<npoints_traj;i++) { for (int j=0;j<npoints_traj;j++) { (*traj_points[i]).add_edge(traj_points[j],muon_dir_reco); } }
  //call update_paths on 0th Point in the vector
  //sort points by their score
  //check whether the 0th Point in the vector is the end point.  If not, loop again
  bool bad_path = false;
  while (true) {
    (*traj_points.front()).update_paths();
    std::sort(traj_points.begin(), traj_points.end(), Comparator());
    if ((*traj_points.front()).id == (endpoint_index+1)) { break; }
    if ((*traj_points.front()).exhausted) {
      std::cout << "unable to traverse particle path" << std::endl;
      bad_path = true;
      break;
    }
  }
  //remove all trajectory Points that are not in the path to the end point
  std::vector<std::vector<double>> trajectory_points_final;
  Point* path_pointer = traj_points[0];
  while(true && !bad_path) {
    trajectory_points_final.push_back((*path_pointer).pos);
    path_pointer = (*path_pointer).prior;
    if (path_pointer == NULL) { break; }
  }

  //add muon start and end vertex to list of trajectory points
  trajectory_points_final.insert(trajectory_points_final.begin(), vtx_muon_end_reco);
  trajectory_points_final.insert(trajectory_points_final.end(),   vtx_muon_start_reco);

  //clean up
  for (int i=0;i<npoints_traj;i++) { delete traj_points[i]; }

  return std::make_tuple(bad_path, trajectory_points_final);
}


//take in point cloud, muon start, and specified segment length
//Forms segments and fits each
std::tuple< std::vector<Track>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> > MCS::form_segs (std::vector<std::vector<double>> vec_points, std::vector<double> muon_start, std::vector<double> muon_end, double seg_length) {
  std::vector<Track> seg_vec;            //vector of each segment of trach
  std::vector<double> distance_vec;      //distance from muon start to middle of each segment
  std::vector<double> angle_vec;         //angle between each segment
  std::vector<double> displacement_vec;  //rms displacement of points within each segment

  //fill track
  Track track;
  int npoints = vec_points.size();
  for(int i=0; i<npoints; i++){ track.add_point(vec_points[i]); }

  //get PCA axes for track
  std::vector<double> com         = {0,0,0};
  std::vector<double> eigenvalues = {0,0,0};
  std::vector<std::vector<double>> axes = MCS::fitPCA(track, com, eigenvalues);
  std::vector<double> aAxis = axes[0];
  std::vector<double> bAxis = axes[1];
  std::vector<double> cAxis = axes[2];
  //flip axis to be along muon direction
  std::vector<double> vec_muon = { muon_end[0]-muon_start[0], muon_end[1]-muon_start[1], muon_end[2]-muon_start[2] };
  if(MCSHelper::dot(vec_muon,aAxis) < 0){ aAxis = MCSHelper::scale(aAxis,-1.); }
  ComparePCAProjection comparator = ComparePCAProjection();
  comparator.sort_points(track.points, aAxis);

  //fill temporary track with remaining points, this gets reduced as they are added to segments
  Track track_remainder = MCS::get_seg(track, aAxis, 1e6);

  std::vector<std::vector<std::vector<double>>> segPCAFit_vec;  //pca fit vectors A,B,C for each 14 cm segment
  std::vector<std::vector<double>> segs_aAxis_vec;  //center of mass position for each segment
  std::vector<std::vector<double>> segCOM_vec;  //center of mass position for each segment
  std::vector<double> angleProjB_vec, angleProjC_vec;

  double currentFirstPointProj = MCSHelper::dot(track.points.front(),aAxis);
  double currentLastPointProj  = MCSHelper::dot(track.points.front(),aAxis);
  double end_proj              = MCSHelper::dot(track.points.back(), aAxis);

  //subdivide into segments of specified length
  //use flexible length to handle gaps
  int iseg = 0;
  while((currentLastPointProj+0.5*seg_length) < end_proj && track_remainder.N>=10){

    segPCAFit_vec.push_back({{0,0,0}, {0,0,0}, {0,0,0}});
    segCOM_vec.push_back(   {0,0,0});
    seg_vec.push_back(Track());

    //Form and fit 14 cm seg with pcafit, update currentFirstPointProj, currentLastPointProj
    bool can_fit = MCS::fitSegPCA(track_remainder, aAxis, seg_vec.back(), segPCAFit_vec.back(), segCOM_vec.back(), seg_length, currentFirstPointProj, currentLastPointProj);
    if (!can_fit) { break; } 
    segs_aAxis_vec.push_back(segPCAFit_vec.back()[0]);
    distance_vec.push_back((currentFirstPointProj+currentLastPointProj)/2 - MCSHelper::dot(muon_start,aAxis));

    //get 3D angle between tracks and 2D projection angles
    if (iseg==0) { angle_vec.push_back(-1); angleProjB_vec.push_back(-1); angleProjC_vec.push_back(-1); }
    else         { MCS::setSegAngles(segPCAFit_vec[iseg-1], segs_aAxis_vec.back(), angle_vec, angleProjB_vec, angleProjC_vec); }

    iseg++;
  }

  return { seg_vec, axes, segs_aAxis_vec, segCOM_vec, distance_vec, angle_vec, angleProjB_vec, angleProjC_vec };
}

double MCS::beta(double gamma) {return sqrt(1.-(1./pow(gamma,2)));}

double MCS::gamma(double KE, double mass) {return KE/mass+1;}	//KE, mass in MeV

double MCS::sigmoid(double x) { return 1. / (1. + std::exp(-x)); }

//This graph gives a function from muon residual range values to total kinetic energy values.
//The data is taken from Atomic Data and Nuclear Data Tables 78: MUON STOPPING POWER AND RANGE TABLES 10 MeV–100 TeV
// http://pdg.lbl.gov/2018/AtomicNuclearProperties/adndt.pdf
// http://pdg.lbl.gov/2017/AtomicNuclearProperties/HTML/liquid_argon.html
// http://pdg.lbl.gov/2017/AtomicNuclearProperties/MUE/muE_liquid_argon.pdf
void MCS::setUKEfromRR() {
  const int un = 20;																	//number of entries in muon KE vs RR table
  //double uCSDA[un]   = { .9937, 1.795, 3.329, 6.605, 10.58, 30.84, 42.50, 67.32, 106.3, 172.5, 238.4, 493.4, 616.3, 855.2, 1202., 1758., 2297., 4359.,  5354.,  7298. };	
  double uCSDA[un]     = { .9833, 1.786, 3.321, 6.598, 10.58, 30.84, 42.50, 67.32, 106.3, 172.5, 238.5, 493.4, 616.3, 855.2, 1202., 1758., 2297., 4359.,  5354.,  7298. };	
  double uEnergy[un]   = {  10.0,  14.0,  20.0,  30.0,  40.0,  80.0, 100.0, 140.0, 200.0, 300.0, 400.0, 800.0, 1000., 1400., 2000., 3000., 4000., 8000., 10000., 14000. };							//Total particle energy in MeV
  double uResRange[un];
  for (int i=0;i<un;i++) { uResRange[i] = uCSDA[i]/rho; }	//cm
  uKEfromRR = new TGraph(un,uResRange,uEnergy);
  uRRfromKE = new TGraph(un,uEnergy,uResRange);
}

//helper function that increases the energy in steps from 10 MeV to 14 GeV 
void MCS::increment_energy (double &e) {
  e *= 1.2;
  e += 2;
}

//helper function that decreases the remaining distance of a muon track
void MCS::decrement_dist (double &x, double xmax) {
  double x_offset = xmax-x;
  x_offset *= 1.1;
  x_offset += 0.5;
  x = xmax - x_offset;
}


//highland formula
double MCS::sigmaH (double T){ return 13.6*(T+Mmu)/T/(T+2*Mmu); }

//helper function to weight the area of each Gaussian using a sigmoid distribution
double MCS::sigmoid(double x, std::vector<double> par) { return par[0] + (par[1]-par[0])*(1 - 1./(1+std::exp(-par[2]*(x/1000.-par[3])))); }

//helper function to modify highland formula based on parameter values
//first 5 parameters define a quartic with the last acting as a scale parameter
double MCS::quartic_decay(double x, std::vector<double> par) {
  double u   = x/par.back()/1000.;
  double val = 0;
  for (int i=0,n=par.size();i<n-1;i++) { val += par[i]*std::pow(u,i); }
  return 1 + val*std::exp(-u);
}

//function to calculate theta_xz PDF parameters
std::vector<double> MCS::pred_theta_xz_pars(double T) {
  double sigma1     = std::sqrt(std::pow(sigmaH(T) * quartic_decay(T,par_sigma1_xz), 2) + std::pow(res_sigma1_xz, 2));
  double sigma2     = std::sqrt(std::pow(sigmaH(T) * quartic_decay(T,par_sigma2_xz), 2) + std::pow(res_sigma2_xz, 2));
  double area_ratio = sigmoid(T, par_ratio_xz);
  return {sigma1, sigma2, area_ratio};
}

//function to calculate theta_yz PDF parameters
std::vector<double> MCS::pred_theta_yz_pars(double T, int vx_index) {
  double sigma1 = std::sqrt(std::pow(sigmaH(T) * quartic_decay(T,par_sigma1_yz[vx_index]), 2) + std::pow(res_sigma1_yz[vx_index], 2));
  double sigma2 = std::sqrt(std::pow(sigmaH(T) * quartic_decay(T,par_sigma2_yz[vx_index]), 2) + std::pow(res_sigma2_yz[vx_index], 2));
  double area_ratio = sigmoid(T, par_ratio_yz[vx_index]);
  return {sigma1, sigma2, area_ratio};
}

//Two-Gaussian PDF that computes the likelihood of an angle given input parameter values
double MCS::double_gaussian(double angle, std::vector<double> pars) {
  double sigma1 = pars[0];
  double sigma2 = pars[1];
  double ratio  = pars[2];
  double gaussian1 = (1./(std::sqrt(2*M_PI)*sigma1)) * std::exp(-0.5*std::pow(angle/sigma1,2));
  double gaussian2 = (1./(std::sqrt(2*M_PI)*sigma2)) * std::exp(-0.5*std::pow(angle/sigma2,2));
  return ratio*gaussian1 + (1-ratio)*gaussian2;
}

//function to computethe likelihood of a given theta_xz angle measurement given a kinetic energy estimate T
double MCS::lnlikelihood_theta_xz(double angle, double T) {
  std::vector<double> pars = pred_theta_xz_pars(T);
  return -std::log(double_gaussian(angle, pars));
}


//function to computethe likelihood of a given theta_yz angle measurement given a kinetic energy estimate T
double MCS::lnlikelihood_theta_yz(double angle, double T, double vx) {

  // Define vx edges and emu edges
  std::vector<double> vx_edges  = {0, 0.1, 0.2, 0.35, 0.75, 1};
  std::vector<double> emu_edges = {600, 950, 1300};

  // Determine the vx_index based on vx
  double vx_abs = std::abs(vx);
  int ivx = 0*(vx_abs>=vx_edges[0] && vx_abs<vx_edges[1]) + 1*(vx_abs>=vx_edges[1] && vx_abs<vx_edges[2]) + 2*(vx_abs>=vx_edges[2] && vx_abs<vx_edges[3]) + 3*(vx_abs>=vx_edges[3] && vx_abs<vx_edges[4]) + 4*(vx_abs>=vx_edges[4] && vx_abs<vx_edges[5]);

  //probability using PDFs from each vx slice
  std::vector<double> pvx = { double_gaussian(angle,pred_theta_yz_pars(T,0)), double_gaussian(angle,pred_theta_yz_pars(T,1)), double_gaussian(angle,pred_theta_yz_pars(T,2)), double_gaussian(angle,pred_theta_yz_pars(T,3)), double_gaussian(angle,pred_theta_yz_pars(T,4)) };

  //apply a different PDF based on vx slice. For highest vx slices, drop down to lower vx slices at large energies
  double probability = 0.0;
  double width = 50;
  double scale1 = sigmoid((T-emu_edges[0])/width);
  double scale2 = sigmoid((T-emu_edges[2])/width);
  if      (ivx == 4) {
    if      (T <  emu_edges[1]) { probability = (1-scale1)*pvx[ivx]   + scale1*pvx[ivx-1]; }
    else if (T >= emu_edges[1]) { probability = (1-scale2)*pvx[ivx-1] + scale2*pvx[ivx-2]; }
  } else if (ivx == 3) {
                                { probability = (1-scale2)*pvx[ivx]   + scale2*pvx[ivx-1]; }
  } else if (ivx <= 2)          { probability = pvx[ivx]; }

  return -std::log(probability);
}

//function to computethe likelihood of a given series of theta_xz and theta_yz angle measurement given an energy estimate E
double MCS::lnlikelihood_track(double* KE, double* par){
  int nsegs = par[0];
  double lnlikelihood = 0;
  for (int i=2; i<nsegs+1; i++) {
    double theta_xz = par[i + nsegs];     // angle in x-projection
    double theta_yz = par[i + 2*nsegs];   // angle in y-projection
    double vx       = par[i + 3*nsegs];   // vx slice
    
    double distance1 = par[i-1];
    double distance2 = par[i];
    double rrtot_guess = uRRfromKE->Eval(KE[0]);
    double rrguess1 = std::max(rrtot_guess - distance1, 1.); //Ensure at least 1 cm of distance
    double rrguess2 = std::max(rrtot_guess - distance2, 1.);
    double keguess1 = uKEfromRR->Eval(rrguess1);
    double keguess2 = uKEfromRR->Eval(rrguess2);
    //double keguess1  = std::max(1.0, uKEfromEX->Interpolate(KE[0], distance1)); // Ensure energy stays above 1 MeV
    //double keguess2  = std::max(1.0, uKEfromEX->Interpolate(KE[0], distance2));
    double keguess   = (keguess1+keguess2)/2 ; // angles is matched to avg energy of the two segments 
	
    //compute likelihood
    double lnl_xz = lnlikelihood_theta_xz(theta_xz, keguess);
    double lnl_yz = lnlikelihood_theta_yz(theta_yz, keguess, vx);
    lnlikelihood += lnl_xz + lnl_yz;
  }
    return lnlikelihood;
}

//Estimate the most likely energy estimate by combining the likelihood predictions from theta_xz and theta_yz
std::vector<double> MCS::estimate_energy(std::vector<double> segs_distance, std::vector<double> segs_angle_x, std::vector<double> segs_angle_y,  std::vector<double> vx_comps){
  // Combine angle and distance vectors for intput into TF1
  segs_distance.insert(segs_distance.begin(), segs_distance.size());  // Add counter as the first entry
  segs_distance.insert(segs_distance.end(), segs_angle_x.begin(), segs_angle_x.end());  // Combine angle_x
  segs_distance.insert(segs_distance.end(), segs_angle_y.begin(), segs_angle_y.end());  // Combine angle_y
  segs_distance.insert(segs_distance.end(), vx_comps.begin(), vx_comps.end());  // Combine vx components

  double emin = 0;
  double emax = 4e3; //4 GeV max estimate

  // Create a TF1 object to minimize the likelihood
  //TF1* f_lnlikelihood2 = new TF1("Negative LnLikelihood of Given Track using Angle-Based MCS", lnlikelihood_track, emin, emax, segs_distance.size());
  TF1* f_lnlikelihood2 = new TF1("Negative LnLikelihood of Given Track using Angle-Based MCS", [this](double* KE, double* par){ return lnlikelihood_track(KE,par); }, emin, emax, segs_distance.size());
  f_lnlikelihood2->SetParameters(&segs_distance[0]);

  // Find the energy that minimizes the likelihood
  double keguess = f_lnlikelihood2->GetMinimumX(emin + 1e-3, emax - 1e-3);

  //define lower and upper bounds
  double keguess_lower  = f_lnlikelihood2->GetMinimumX(                      emin + 1e-3,  keguess*0.8);
  double keguess_higher = f_lnlikelihood2->GetMinimumX(std::min(keguess*1.2, emax - 2e-3), emax - 1e-3);

  //get likelihood at e_guess, e_guess_lower, e_guess_higher
  double l_keguess        = std::exp(-lnlikelihood_track(&keguess,       &segs_distance[0]));
  double l_keguess_lower  = std::exp(-lnlikelihood_track(&keguess_lower, &segs_distance[0]));
  double l_keguess_higher = std::exp(-lnlikelihood_track(&keguess_higher,&segs_distance[0]));

  //copmute the ambiguity score as the highest ratio
  double ambiguity_score = std::max(l_keguess_lower/l_keguess ,l_keguess_higher/l_keguess );

  // Clean up
  delete f_lnlikelihood2;
  return {keguess, ambiguity_score};
} 
