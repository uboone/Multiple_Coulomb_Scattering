#include "TROOT.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "src/mcs.h"
#include "src/mcs.cxx"

int main (int argc, char *argv[]) {

  //Read in data
  std::string filename = "data/simulated_event.root";
  TFile* file = TFile::Open(filename.c_str());
  TMatrixD* traj_points_mat = (TMatrixD*)file->Get("spacepoints_evt0");
  TVectorD* start_vtx = (TVectorD*)file->Get("reco_start_evt0");
  TVectorD* end_vtx   = (TVectorD*)file->Get("reco_end_evt0");
  std::vector<double> vtx_muon_start_reco = { (*start_vtx)[0], (*start_vtx)[1], (*start_vtx)[2] };
  std::vector<double> vtx_muon_end_reco   = {   (*end_vtx)[0],   (*end_vtx)[1],   (*end_vtx)[2] };
  int npoints_trajectory = traj_points_mat->GetNrows();
  std::vector<std::vector<double>> trajectory_points_initial;
  for (int i=0;i<npoints_trajectory;i++) { trajectory_points_initial.push_back({ (*traj_points_mat)(i,0), (*traj_points_mat)(i,1), (*traj_points_mat)(i,2) }); }

  //Create MCS class
  MCS mcs;

  //run over event(s)
  //In this example there is one simulated muon with true energy 0.735 GeV, which MCS reconstructs at 0.699 GeV
  //The algortihm uses three inputs: the reconstructed start and end muon vertices, and a vector of 3D spacepoints following the muon path
  //Any 3D pointcloud should suffice, but it will help to remove artifacts such as delta rays and crossing tracks
  mcs.run(vtx_muon_start_reco, vtx_muon_end_reco, trajectory_points_initial);

  //clean up memory
  mcs.cleanUp();
}
