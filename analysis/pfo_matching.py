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
parser.add_argument('-i-', '--inputFile', help='--inputFile output_reco.slcio', 
                  type=str, default='output_reco.slcio')

# Output file
parser.add_argument('-o', '--outputFile', help='--outputFile pfo_matching.root', 
                  type=str, default='pfo_matching.root')
args = parser.parse_args()

hists = []

# Initialize histograms
fMatchedPt = TH1F('matched_pt', 'Matched PFO True Pt', 30, 0, 300)
hists.append(fMatchedPt)
fMatchedEta = TH1F('matched_eta', 'Matched PFO True Eta', 25, -3, 3)
hists.append(fMatchedEta)
fAllPt = TH1F('all_pt', 'All True Pt', 30, 0, 300)
hists.append(fAllPt)
fAllEta = TH1F('all_eta', 'All True Eta', 25, -3, 3)
hists.append(fAllEta)

for hist in hists:
  hist.SetDirectory(0)
  
to_process = []

if os.path.isdir(args.inputFile):
  for r, d, f in os.walk(args.inputFile):
    for file in f:
      to_process.append(os.path.join(r, file))
else:
  to_process.append(args.inputFile)

def eta(theta):
  return (-math.log(math.tan(theta/2)))

# Open input file(s)
for file in to_process:
  reader = IOIMPL.LCFactory.getInstance().createLCReader()
  reader.open(file)

  # Loop through Events
  for ievt, event in enumerate(reader):

    # Get collections
    pfos = event.getCollection('PandoraPFOs')
    mcs = event.getCollection('MCParticle')

    # Loop over MCParticle collection
    for mc in mcs:
      
      # Select charged pions
      if (abs(mc.getPDG()) != 211):
        continue

      # Get mc pion kinematic variables
      mcE = mc.getEnergy()
      mcMom = mc.getMomentum()
      mcPt = math.sqrt(mcMom[0]**2+mcMom[1]**2)

      # Cut on pT?

      mcTheta = math.acos(mcMom[2]/(math.sqrt(mcPt**2+mcMom**2)))
      mcEta = eta(mcTheta)

      # Cut on theta?

      mcPhi = math.acos(mcMom[0]/mcPt)

      # Initialize large minimum dR
      min_dR = 1000

      # Loop through pfo collection
      for pfo in pfos:

        # Select pions
        if (abs(pfo.getType() != 211):
          continue

        # Get cluster collection
        clusters = pfo.getClusters()

        # Initialize zero minimum dE

        # Loop over clusters
        for cluster in clusters:
          # Calculate energy difference
          E_diff = cluster.getEnergy() - mcE
          if (math.fabs(E_diff) > dE):
            dE = E_diff

        # Get pfo pion kinematic variables
        pfoMom = pfo.getMomentum()
        pfoPt = math.sqrt(pfoMom[0]**2+pfoMom[1]**2)
        pfoTheta = math.acos(pfoMom[2]/(math.sqrt(pfoPt**2+pfoMom**2)))
        pfoEta = eta(pfoTheta)
        pfoPhi = math.acos(pfoMom[0]/pfoPt)

        # Calculate dR
        dR = math.sqrt((mcPhi-pfoPhi)**2+(mcEta-pfoEta)**2)
        if (dR < min_dR):
          min_dR = dR

      if (min_dR < 0.25):
        fMatchedPt.Fill(mcPt)
        fMatchedEta.Fill(mcEta)

      fAllPt.Fill(mcPt)
      fAllEta.Fill(mcEta)
      
  reader.close()

# Create charged pion pt efficiency hist
fPiPtEff = fMatchedPt.Clone('pt_eff')
fPiPtEff.Divide(fPiPtEff, fAllPt, 1, 1, 'B')
fPiPtEff.SetLineColor(6)
fPiPtEff.SetLineWidth(2)
fPiPtEff.SetTitle('Charged Pion Efficiency vs Pt')
fPiPtEff.GetXaxis().SetTitle('pT [GeV/c]')
fPiPtEff.GetYaxis().SetTitle("#epsilon")
fPiPtEff.SetStats(0)
hists.append(fPiPtEff)

# Create charged pion eta efficiency hist
fPiEtaEff = fMatchedEta.Clone('eta_eff')
fPiEtaEff.Divide(fPiEtaEff, fAllEta, 1, 1, 'B')
fPiEtaEff.SetLineColor(7)
fPiEtaEff.SetLineWidth(2)
fPiEtaEff.SetTitle('Charged Pion Efficiency vs Eta')
fPiEtaEff.GetXaxis().SetTitle('#eta')
fPiEtaEff.GetYaxis().SetTitle('#epsilon')
fPiEtaEff.SetStats(0)
hists.append(fPiEtaEff)

output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
  hist.Write()
output_file.Close()

for hist in hists:
  filename = hist.Name() + '.png'
  canvas = TCanvas()
  hist.Draw()
  canvas.SaveAs(filename)
