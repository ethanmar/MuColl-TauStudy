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
                  type=str, default='output_digi_light.slcio')

# Output file
parser.add_argument('-o', '--outputFile', help='--outputFile pi_ana.root', 
                  type=str, default='pi_ana.root')
args = parser.parse_args()

hists = []

# Initialize histograms
fPDG = TH1F('pdg', 'MC Particle PDG', 1250, 0, 2500)
hists.append(fPDG)
fType = TH1F('type', 'Reco PFO Type', 1250, 0, 2500)
hists.append(fType)
fNMCPions = TH1F('n_mc_pions', 'Number of MC Pions', 400, 0, 400)
hists.append(fNMCPions)
fNRecoPions = TH1F('n_reco_pions', 'Number of Reco Pions', 400, 0, 400)
hists.append(fNRecoPions)
fMCPiPt = TH1F('mc_pi_pt', 'MC Pion Pt', 20, 0, 110)
fMCPiPt.SetXTitle('Pt [GeV/c]')
hists.append(fMCPiPt)
fMCPi0Pt = TH1F('mc_pi0_pt', 'MC Neutral Pion Pt', 20, 0, 110)
fMCPi0Pt.SetXTitle('Pt [GeV/c]')
hists.append(fMCPi0Pt)
fRecoPiPt = TH1F('reco_pi_pt', 'Reco Pion Pt', 20, 0, 110)
fRecoPiPt.SetXTitle('Pt [GeV/c]')
hists.append(fRecoPiPt)
fRecoPi0Pt = TH1F('reco_pi0_pt', 'Reco Neutral Pion Pt', 20, 0, 110)
fRecoPi0Pt.SetXTitle('Pt [GeV/c]')
hists.append(fRecoPi0Pt)

for hist in hists:
  hist.SetDirectory(0)

pfoTypes = {}
  
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

  # Loop through Events
  for ievt, event in enumerate(reader):
    # print('Entered event')
    mc_pis = []
    mc_pi0s = []
    reco_pis = []
    reco_pi0s = []
    
    mc_particles = event.getCollection('MCParticle')

    # Loop through MCParticles
    for mc_particle in mc_particles:
      #print('Entered MCParticle')
      pdg = abs(mc_particle.getPDG())
      parents = mc_particle.getParents()
      if (len(parents) != 0):
        motherID = parents[-1].getPDG()
      else:
        motherID = 0

      # Tag charged and neutral pions
      if (pdg == 211 and motherID == 15):
        mc_pis.append(mc_particle)
      elif (pdg == 111 and motherID == 15):
        mc_pi0s.append(mc_particle)

      # Fill PDG hist
      fPDG.Fill(pdg)

    # Fill NPions hist
    fNMCPions.Fill(len(mc_pis)+len(mc_pi0s))

    # Loop over charged pions
    for mc_pi in mc_pis:
      p = mc_pi.getMomentum()
      px = p[0]
      py = p[1]
      pt = math.sqrt(px**2 + py**2)

      # Fill MC charged pion pt hist
      fMCPiPt.Fill(pt)

    # Loop over neutral pions
    for mc_pi0 in mc_pi0s:
      p = mc_pi0.getMomentum()
      px = p[0]
      py = p[1]
      pt = math.sqrt(px**2 + py**2)

      # Fill MC neutral pion pt hist
      fMCPi0Pt.Fill(pt)

    pfos = event.getCollection('PandoraPFOs')

    # Loop over PFOs (reconstructed particles)
    for pfo in pfos:
      # print('Entered PFO')
      _type = abs(pfo.getType())

      if str(_type) in pfoTypes:
        pfoTypes[str(_type)] += 1
      else:
        pfoTypes[str(_type)] = 1
      
      # Tag charged and neutral pions
      if _type == 211:
        reco_pis.append(pfo)
      elif _type == 111:
        reco_pi0s.append(pfo)

      # Fill PFO type hist (reconstructed particle type)  
      fType.Fill(_type)

    # Fill NPions hist  
    fNRecoPions.Fill(len(reco_pis) + len(reco_pi0s))

    # Loop over charged pions
    for reco_pi in reco_pis:
      p = reco_pi.getMomentum()
      px = p[0]
      py = p[1]
      pt = math.sqrt(px**2 + py**2)

      # Fill reco charged pion pt hist
      fRecoPiPt.Fill(pt)

    # Loop over neutral pions
    for reco_pi0 in reco_pi0s:
      p = reco_pi0.getMomentum()
      px = p[0]
      py = p[1]
      pt = math.sqrt(px**2 + py**2)

      # Fill reco neutral pion pt hist
      fRecoPi0Pt.Fill(pt)  
      
  reader.close()

# Create charged pion pt efficiency hist
fPiPtEff = fRecoPiPt.Clone('pi_eff')
fPiPtEff.Divide(fPiPtEff, fMCPiPt, 1, 1, 'B')
fPiPtEff.SetLineColor(6)
fPiPtEff.SetLineWidth(2)
fPiPtEff.SetTitle('Charged Pion Efficiency vs Pt')
fPiPtEff.GetXaxis().SetTitle('pT [GeV/c]')
fPiPtEff.GetYaxis().SetTitle("#epsilon")
hists.append(fPiPtEff)

# Create neutral pion pt efficiency hist
fPi0PtEff = fRecoPi0Pt.Clone('pi0_eff')
fPi0PtEff.Divide(fPi0PtEff, fMCPi0Pt, 1, 1, 'B')
fPi0PtEff.SetLineColor(7)
fPi0PtEff.SetLineWidth(2)
fPi0PtEff.SetTitle('Neutral Pion Efficiency vs Pt')
fPi0PtEff.GetXaxis().SetTitle('pT [GeV/c]')
fPi0PtEff.GetYaxis().SetTitle('#epsilon')
hists.append(fPi0PtEff)

output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
  hist.Write()
output_file.Close()


for hist in hists:
  filename = hist.GetTitle() + '.png'
  canvas = TCanvas()
  hist.Draw()
  canvas.SaveAs(filename)

with open('pfo_types.txt', 'w') as file:
  file.write('PDG/Type of Reconstructed PFOs')
  for key, value in pfoTypes.items():
    print(f'PDG: {key}, Total Number: {value}', file=file)
