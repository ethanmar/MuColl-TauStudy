from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array 
import os
import fnmatch

# Command line arguments
parser = ArgumentParser()

# Input file 
parser.add_argument('-i', '--inputFile', help='--inputFile output_reco.slcio', 
                    type=str, default='output_reco.slcio')

# Output file
parser.add_argument('-o', '--outputFile', help='--outputFile cluster_ana.root',
                    type=str, default='cluster_ana.root')

args = parser.parse_args()

hists = []

# Initialize histograms
fEnergy = TH1F('energy', 'Energy', 100, 0, 100)
hists.append(fEnergy)
fPhi = TH1F('phi', 'Phi', 180, 0, 360)
hists.append(fPhi)
fTheta = TH1F('theta', 'Theta', 180, 0, 360)
hists.append(fTheta)
fNHits = TH1F('nHits', 'Number of Calo Hits', 50, 0, 50)
hists.append(fNHits)

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
    
    # Initialize pandora cluster collection
    cluster_coll = event.getCollection('PandoraClusters')

    # Loop through pandora cluster collection
    for cluster in cluster_coll:

      energy = cluster.getEnergy()
      phi = cluster.getIPhi()
      theta = cluster.getITheta()
      caloHits = cluster.getCalorimeterHits()
      nHits = len(caloHits)

      fEnergy.Fill(energy)
      fPhi.Fill(phi)
      fTheta.Fill(theta)
      fNHits.Fill(nHits)

  # Close input file
  reader.close()

# Write to output file
output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
  hist.Write()
output_file.Close()

# Plot histograms
for hist in hists:
  filename = hist.GetName() + '.png'
  canvas = TCanvas()
  hist.Draw()
  canvas.SaveAs(filename)
