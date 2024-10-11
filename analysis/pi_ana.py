from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
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

hists = []

hMC_Pdg = TH1F('MC_Pdg', 'MC_Pdg', 3000, 0, 3000)
hists.append(hMC_Pdg)
hPFO_Type = TH1F('PFO_type', 'PFO_type', 3000, 0, 3000)
hists.append(hPFO_Type)
hSel_PFO_Type = TH1F('Sel_PFO_type', 'Sel_PFO_type', 3000, 0, 3000)
hists.append(hSel_PFO_Type)
hN_MCPions = TH1F('N_MCPions', 'N_MCPions', 400, 0, 400)
hists.append(hN_MCPions)
hN_Pions = TH1F('N_Pions', 'N_Pions', 400, 0, 400)
hists.append(hN_Pions)
hN_Sel_Pions = TH1F('N_Sel_Pions', 'N_Sel_Pions', 400, 0, 400)
hists.append(hN_Sel_Pions)
hMC_Pt = TH1F('MC_Pt', 'MC_Pt', 1000, 0, 250)
hists.append(hMC_Pt)
hPFO_Pt = TH1F('PFO_Pt', 'PFO_Pt', 1000, 0, 250)
hists.append(hPFO_Pt)
hSel_PFO_Pt = TH1F('Sel_PFO_Pt', 'Sel_PFO_Pt', 1000, 0, 250)
hists.append(hSel_PFO_Pt)


for hist in hists:
  hist.SetDirectory(0)

to_process = []

if os.path.isdir(args.inputFile):
  for r, d, f in os.walk(args.inputFile):
    for file in f:
      to_process.append(os.path.join(r, file))
else:
  to_process.append(args.inputFile)

  
for file in to_process:
  reader = IOIMPL.LCFactory.getInstance().createLCReader()
  reader.open(file)
  for ievt, event in enumerate(reader):
    mc_pions = []
    pfo_pions = []
    sel_pfo_pions = []
    
    mc_particles = event.getCollection('MCParticle')
    for mc_particle in mc_particles:
      mc_pdg = abs(mc_particle.getPDG())

      if mc_pdg == 211:
        mc_pions.append(mc_particle)

      hMC_Pdg.Fill(mc_pdg)

    hN_MCPions.Fill(len(mc_pions))
      
    for mc_pion in mc_pions:
      mc_p = mc_pion.getMomentum()
      mc_px = mc_p[0]
      mc_py = mc_p[1]
      mc_pt = math.sqrt(mc_px**2 + mc_py**2)

      hMC_Pt.Fill(mc_pt)

    pfos = event.getCollection('PandoraPFOs')
    for pfo in pfos:
      pfo_type = abs(pfo.getType())

      if pfo_type == 211:
        pfo_pions.append(pfo)

      hPFO_Type.Fill(pfo_type)

    hN_Pions.Fill(len(pfo_pions))

    for pfo_pion in pfo_pions:
      pfo_p = pfo_pion.getMomentum()
      pfo_px = pfo_p[0]
      pfo_py = pfo_p[1]
      pfo_pt = math.sqrt(pfo_px**2 + pfo_py**2)

      hPFO_Pt.Fill(pfo_pt)

    sel_pfos = event.getCollection('SelectedPandoraPFOs')
    for sel_pfo in sel_pfos:
      sel_pfo_type = abs(sel_pfo.getType())

      if sel_pfo_type == 211:
        sel_pfo_pions.append(sel_pfo)

      hSel_PFO_Type.Fill(sel_pfo_type)

    hN_Sel_Pions.Fill(len(sel_pfo_pions))

    for sel_pfo_pion in sel_pfo_pions:
      sel_pfo_p = sel_pfo_pion.getMomentum()
      sel_pfo_px = sel_pfo_p[0]
      sel_pfo_py = sel_pfo_p[1]
      sel_pfo_pt = math.sqrt(sel_pfo_px**2 + sel_pfo_py**2)

      hSel_PFO_Pt.Fill(sel_pfo_pt)

      
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


  
