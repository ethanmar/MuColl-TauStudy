from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile
from math import *
from argparse import ArgumentParser
from array import array
import os
import fnmatch

parser = ArgumentParser()

parser.add_argument('-i-', '--inputFile', help='--inputFile output_digi_light.slcio', 
                  type=str, default='output_digi_light.slcio')
parser.add_argument('-o', '--outputFile', help='--outputFile pdg_ana.root', 
                  type=str, default='pdg_ana.root')
args = parser.parse_args()

h_MC_Pdg = TH1F('MC_Pdg', 'MC_Pdg', 3000, 0, 3000)

h_MC_Pdg.SetDirectory(0)

to_process = []

if os.path.isdir(args.inputFile):
  for r, d, f in os.walk(args.inputFile):
    for file in f:
      to_process.append(os.path.join(r, file))
else:
  to_process.append(args.inputFile)

filenum = 0

for file in to_process:
  reader = IOIMPL.LCFactory.getInstance().createLCReader()
  reader.open(file)
  filenum = filenum + 1
  for ievt, event in enumerate(reader):
    mcpCollection = event.getCollection('MCParticle')
    particle = mcpCollection.getElementAt(i)
    pdg = particle.getPDG()
    h_MC_Pdg.Fill(abs(pdg))
  reader.close()

output_file = TFile(args.outputFile, 'RECREATE')
h_MC_Pdg.Write()
output_file.Close()
  
