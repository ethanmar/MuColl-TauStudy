from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

args = parser.parse_args()

# Initialize histograms

hists = []

hNPFOs = TH1F('nPFOs', 'Number of Reconstructed PFOs', 35, 0, 35)
hists.append(hNPFOs)
hEnergyPFOs_11 = TH1F('energyPFOs_pdg11', 'Energy of Reconstructed PFOs (PDG 11)', 50, 0, 10)
hists.append(hEnergyPFOs_11)
hPtPFOs_11 = TH1F('ptPFOs_pdg11', 'Transverse Momentum of Reconstructed PFOs (PDG 11)', 50, 0, 10)
hists.append(hPtPFOs_11)
hEnergyRatioPFOs_11 = TH1F('energyRatioPFOs_pdg11', 'Energy Ratio of Reconstructed PFOs (PDG 11)', 50, 0, 1)
hists.append(hEnergyRatioPFOs_11)
hEnergyPFOs_13 = TH1F('energyPFOs_pdg13', 'Energy of Reconstructed PFOs (PDG 13)', 50, 0, 10)
hists.append(hEnergyPFOs_13)
hPtPFOs_13 = TH1F('ptPFOs_pdg13', 'Transverse Momentum of Reconstructed PFOs (PDG 13)', 50, 0, 10)
hists.append(hPtPFOs_13)
hEnergyRatioPFOs_13 = TH1F('energyRatioPFOs_pdg13', 'Energy Ratio of Reconstructed PFOs (PDG 13)', 50, 0, 1)
hists.append(hEnergyRatioPFOs_13)
hEnergyPFOs_22 = TH1F('energyPFOs_pdg22', 'Energy of Reconstructed PFOs (PDG 22)', 50, 0, 10)
hists.append(hEnergyPFOs_22)
hPtPFOs_22 = TH1F('ptPFOs_pdg22', 'Transverse Momentum of Reconstructed PFOs (PDG 22)', 50, 0, 10)
hists.append(hPtPFOs_22)
hEnergyRatioPFOs_22 = TH1F('energyRatioPFOs_pdg22', 'Energy Ratio of Reconstructed PFOs (PDG 22)', 50, 0, 1)
hists.append(hEnergyRatioPFOs_22)
hEnergyPFOs_211 = TH1F('energyPFOs_pdg211', 'Energy of Reconstructed PFOs (PDG 211)', 50, 0, 10)
hists.append(hEnergyPFOs_211)
hPtPFOs_211 = TH1F('ptPFOs_pdg211', 'Transverse Momentum of Reconstructed PFOs (PDG 211)', 50, 0, 10)
hists.append(hPtPFOs_211)
hEnergyRatioPFOs_211 = TH1F('energyRatioPFOs_pdg211', 'Energy Ratio of Reconstructed PFOs (PDG 211)', 50, 0, 1)
hists.append(hEnergyRatioPFOs_211)
hEnergyPFOs_2112 = TH1F('energyPFOs_pdg2112', 'Energy of Reconstructed PFOs (PDG 2112)', 50, 0, 10)
hists.append(hEnergyPFOs_2112)
hPtPFOs_2112 = TH1F('ptPFOs_pdg2112', 'Transverse Momentum of Reconstructed PFOs (PDG 2112)', 50, 0, 10)
hists.append(hPtPFOs_2112)
hEnergyRatioPFOs_2112 = TH1F('energyRatioPFOs_pdg2112', 'Energy Ratio of Reconstructed PFOs (PDG 2112)', 50, 0, 1)
hists.append(hEnergyRatioPFOs_2112)
    
# Detach histograms from file/directory
for hist in hists:
    hist.SetDirectory(0)

# Check if input file is a directory or a single file    
to_process = []

if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)

n_taus = 0
    
# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        taus = event.getCollection('RecoTaus')
        taus_nparticles = event.getCollection('RecoTausNParticles')
        taus_nchargedtrks = event.getCollection('RecoTausNChargedTrks')
        taus_merge = event.getCollection('RecoTausMerge')
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')

        # Loop through rec taus which passed
        for tau in taus:
            n_taus += 1
            tau_energy = tau.getEnergy()
            pfos = tau.getParticles()
            npfos = len(pfos)
            hNPFOs.Fill(npfos)
            for pfo in pfos:
                pfo_type = abs(pfo.getType())
                E = pfo.getEnergy()
                E_ratio = E/tau_energy
                px = pfo.getMomentum()[0]
                py = pfo.getMomentum()[1]
                pz = pfo.getMomentum()[2]
                pt = math.sqrt(px**2 + py**2)
                p = math.sqrt(pt**2 + pz**2)

                if pfo_type == 11:
                    hEnergyPFOs_11.Fill(E)
                    hPtPFOs_11.Fill(pt)
                    hEnergyRatioPFOs_11.Fill(E_ratio)

                elif pfo_type == 13:
                    hEnergyPFOs_13.Fill(E)
                    hPtPFOs_13.Fill(pt)
                    hEnergyRatioPFOs_13.Fill(E_ratio)
                                    
                elif pfo_type == 22:
                    hEnergyPFOs_22.Fill(E)
                    hPtPFOs_22.Fill(pt)
                    hEnergyRatioPFOs_22.Fill(E_ratio)
                                    
                elif pfo_type == 211:
                    hEnergyPFOs_211.Fill(E)
                    hPtPFOs_211.Fill(pt)
                    hEnergyRatioPFOs_211.Fill(E_ratio)
                                    
                elif pfo_type == 2112:
                    hEnergyPFOs_2112.Fill(E)
                    hPtPFOs_2112.Fill(pt)
                    hEnergyRatioPFOs_2112.Fill(E_ratio)
                                
        # Loop through reco taus which failed number of charged tracks
        for tau_nchargedtrks in taus_nchargedtrks:
            n_taus += 1
            tau_energy = tau_nchargedtrks.getEnergy()
            pfos = tau_nchargedtrks.getParticles()
            npfos = len(pfos)
            hNPFOs.Fill(npfos)
            for pfo in pfos:
                pfo_type = abs(pfo.getType())
                E = pfo.getEnergy()
                E_ratio = E/tau_energy
                px = pfo.getMomentum()[0]
                py = pfo.getMomentum()[1]
                pz = pfo.getMomentum()[2]
                pt = math.sqrt(px**2 + py**2)
                p = math.sqrt(pt**2 + pz**2)

                if pfo_type == 11:
                    hEnergyPFOs_11.Fill(E)
                    hPtPFOs_11.Fill(pt)
                    hEnergyRatioPFOs_11.Fill(E_ratio)

                elif pfo_type == 13:
                    hEnergyPFOs_13.Fill(E)
                    hPtPFOs_13.Fill(pt)
                    hEnergyRatioPFOs_13.Fill(E_ratio)
                                
                elif pfo_type == 22:
                    hEnergyPFOs_22.Fill(E)
                    hPtPFOs_22.Fill(pt)
                    hEnergyRatioPFOs_22.Fill(E_ratio)
                                
                elif pfo_type == 211:
                    hEnergyPFOs_211.Fill(E)
                    hPtPFOs_211.Fill(pt)
                    hEnergyRatioPFOs_211.Fill(E_ratio)
                                
                elif pfo_type == 2112:
                    hEnergyPFOs_2112.Fill(E)
                    hPtPFOs_2112.Fill(pt)
                    hEnergyRatioPFOs_2112.Fill(E_ratio)
                        
        # Loop through reco taus which failed merge
        for tau_merge in taus_merge:
            n_taus += 1
            tau_energy = tau_merge.getEnergy()
            pfos = tau_merge.getParticles()
            npfos = len(pfos)
            hNPFOs.Fill(npfos)
            for pfo in pfos:
                pfo_type = abs(pfo.getType())
                E = pfo.getEnergy()
                E_ratio = E/tau_energy
                px = pfo.getMomentum()[0]
                py = pfo.getMomentum()[1]
                pz = pfo.getMomentum()[2]
                pt = math.sqrt(px**2 + py**2)
                p = math.sqrt(pt**2 + pz**2)
                            
                if pfo_type == 11:
                    hEnergyPFOs_11.Fill(E)
                    hPtPFOs_11.Fill(pt)
                    hEnergyRatioPFOs_11.Fill(E_ratio)

                elif pfo_type == 13:
                    hEnergyPFOs_13.Fill(E)
                    hPtPFOs_13.Fill(pt)
                    hEnergyRatioPFOs_13.Fill(E_ratio)
                                
                elif pfo_type == 22:
                    hEnergyPFOs_22.Fill(E)
                    hPtPFOs_22.Fill(pt)
                    hEnergyRatioPFOs_22.Fill(E_ratio)
                                
                elif pfo_type == 211:
                    hEnergyPFOs_211.Fill(E)
                    hPtPFOs_211.Fill(pt)
                    hEnergyRatioPFOs_211.Fill(E_ratio)
                                
                elif pfo_type == 2112:
                    hEnergyPFOs_2112.Fill(E)
                    hPtPFOs_2112.Fill(pt)
                    hEnergyRatioPFOs_2112.Fill(E_ratio)
                                
        # Loop through reco taus which failed number of particles
        for tau_nparticles in taus_nparticles:
            n_taus += 1
            tau_energy = tau_nparticles.getEnergy()
            pfos = tau_nparticles.getParticles()
            npfos = len(pfos)
            hNPFOs.Fill(npfos)
            for pfo in pfos:
                pfo_type = abs(pfo.getType())
                E = pfo.getEnergy()
                E_ratio = E/tau_energy
                px = pfo.getMomentum()[0]
                py = pfo.getMomentum()[1]
                pz = pfo.getMomentum()[2]
                pt = math.sqrt(px**2 + py**2)
                p = math.sqrt(pt**2 + pz**2)
                            
                if pfo_type == 11:
                    hEnergyPFOs_11.Fill(E)
                    hPtPFOs_11.Fill(pt)
                    hEnergyRatioPFOs_11.Fill(E_ratio)

                elif pfo_type == 13:
                    hEnergyPFOs_13.Fill(E)
                    hPtPFOs_13.Fill(pt)
                    hEnergyRatioPFOs_13.Fill(E_ratio)
                    
                elif pfo_type == 22:
                    hEnergyPFOs_22.Fill(E)
                    hPtPFOs_22.Fill(pt)
                    hEnergyRatioPFOs_22.Fill(E_ratio)
                                
                elif pfo_type == 211:
                    hEnergyPFOs_211.Fill(E)
                    hPtPFOs_211.Fill(pt)
                    hEnergyRatioPFOs_211.Fill(E_ratio)
                                
                elif pfo_type == 2112:
                    hEnergyPFOs_2112.Fill(E)
                    hPtPFOs_2112.Fill(pt)
                    hEnergyRatioPFOs_2112.Fill(E_ratio)
                    
    # Close file
    reader.close()

print(f'Number of taus: {n_taus}')
    
# Write to output file
output_file = TFile('pfo_ana.root', 'RECREATE')
for hist in hists:
    hist.Write()
output_file.Close()

# Draw hists and save as PNG
for hist in hists:
    filename = hist.GetName() + '.png'
    canvas = TCanvas()
    hist.Draw()
    canvas.SaveAs(filename)
