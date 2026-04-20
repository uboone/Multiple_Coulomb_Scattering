#Multiple Coulomb Scattering in MicroBooNE

This example code demonstrates the newly developed Multiple Coulomb Scattering algorithm in <arxiv_link_here>.

Simulation from a single muon in the MicroBooNE detector is provided in data/simulated_event.root, including the reconstructed muon start and end vertices, as well as a vector of 3D points along the muon path. These are used as inputs for the MCS algorithm to demonstrate its performance. The muon was simulated with an energy of 0.735 GeV and is reconstructed with the MCS algorithm to have an estimated energy of 0.699 GeV.

To use the algorithm, first run "make", then run the executable "mcs_example" that is generated.
