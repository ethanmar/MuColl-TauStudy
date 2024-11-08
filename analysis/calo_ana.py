from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array 
import os
import fnmatch

# Command line arguments
parser = ArgumenrParser()

# Input file 
parser.add_argument('-i', '--inputFile', help='--inputFile tau_sim.slcio', 
                    type=str, default='tau_sim.slcio')

# Output file
parser.add_argument('-o', '--outputFile', help='--outputFile cal_ana.root',
                    type=str, default='cal_ana.root')

args = parser.parse_args()

hists = []

# Initialize histograms
fCaloHitPDG = TH1F('calo_hit_pdg', 'calo_hit_pdg', 1250, 0, 2500)
hists.append(fCaloHitPDG)

for hist in hists:
  hist.SetDirectory(0)

to_process = []

if os.path.isdir(args.inputFile):
  for r, d, f in os.walk(args.inputFile):
    for file in f:
      to_process.append(os.path.join(r, file))
else:
  to_process.append(args.inputFile)

# Open input file(s)
for file in to_process:
  reader = IOIMPL.LCFactory.getInstance().createLCReader()
  reader.open(file)

  # Loop through events
  for ievt, event in enumerate(reader):
    # Get SimCalorimeterHit collection
    sim_calo_hits = event.getCollection("SimCalorimeterHits")
    # Loop over simulated calorimeter hits
    for sim_calo_hit in sim_calo_hits:
    
      # Get number of MCParticle contributions to hit
      n_contributions = sim_calo_hit.getNMContributions()

      # Loop over each MCParticle contribution
      for n in range(n_contributions):
        mc_particle = sim_calo_hit.getParticleCont(n)
        pdg = abs(mc_particle.getPDG())

        # Fill calo hit pdg hist
        fCaloHitPDG.Fill(pdg)

  # Close input file
  reader.close()

output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
  hist.Write()
output_file.Close()

for hist in hists:
  filename = hist.GetTitle() + '.png'
  canvas = TCanvas()
  hist.Draw()
  canvas.SaveAs(filename)
