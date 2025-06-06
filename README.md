# Instructions for Complete Workflow of Tau (Decay Product) Event Analysis:
This workflow assumes that you are in the `MuColl-TauStudy` directory:

`cd <your path to>/MuColl-TauStudy`

### Tau Event Generation:
Generate `n` tau events and write them to `gen_output.slcio`:

`python generation/lcio_tau_gun.py --nEvents n --outputFile gen_output.slcio`

This creates one tau MCParticle per event with uniform distributions of 20 <= pT <= 320 GeV/c, 10 <= theta <= 170 deg, and 0 <= phi <= 2pi rad, which can be manually adjusted.

### Event Simulation in MAIA Detector Geometry:
Before simulating the events, change the `MUCOLL_GEO` environment variable to the path of the MAIA detector geometry file:

`export MUCOLL_GEO=/your/path/to/MuColl-TauStudy/detector_geometry/MAIA_v0/MAIA_v0.xml`

Now simulate the decays and material interactions of the generated tau events in the MAIA detector:

`ddsim --steeringFile simulation/steer_sim_tau_gun.py --inputFile gen_output.slcio --outputFile sim_output.slcio`

### Digitization of Simulated Hits:
Convert simulated hits into realistic hits:

`k4run digitisation/k4run/digi_steer.py --LcioEvent.Files sim_output.slcio --outputFile digi_output`

This produces two output files: `digi_output.slcio` and `digi_output_light.slcio`. The former contains all collections produced from the digitization, while the latter contains only a small subset.

### Reconstruction of Tau Decay Products:
Reconstruct tracks, clusters, and particles of tau decay products (plus other simulated particles):

`cp -a reconstruction/k4run/PandoraSettings/ ./`

`k4run reconstruction/k4run/reco_steer.py --LcioEvent.Files digi_output.slcio --MatFile ${ACTS_MatFile} -- TGeoFile ${ACTS_TGeoFile}`

Like the digitization steps, this produces two output files: `reco_output.slcio` and `reco_output_light.slcio`. The former contains all collections produced from the reconstruction, while the latter contains only a small subset.
