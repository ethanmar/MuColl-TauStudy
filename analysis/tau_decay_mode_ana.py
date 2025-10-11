from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt

from tau_mc_link import getDecayMode

# Function to get visible properties of mc tau
def getVisibleProperties(mcTau):

    E_vis = 0
    px_vis = 0
    py_vis = 0
    pz_vis = 0
    
    # Loop through mcDaughters
    mcDaughters = mcTau.getDaughters()

    for mcDaughter in mcDaughters:

        # Ignore neutrinos
        pdg = abs(mcDaughter.getPDG())
        if (pdg ==12 or pdg == 14 or pdg == 16):
            continue

        # Sum visible observables
        E_vis += mcDaughter.getEnergy()
        px_vis += mcDaughter.getMomentum()[0]
        py_vis += mcDaughter.getMomentum()[1]
        pz_vis += mcDaughter.getMomentum()[2]

    vis_properties = [E_vis, px_vis, py_vis, pz_vis]
    return vis_properties

# Script to analyze individual tau decay modes

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

args = parser.parse_args()

# Loop over decay modes

decay_modes = ['1P0N', '3P0N']
decay_mode_nums = [0, 4]

for n, decay_mode in enumerate(decay_modes):

    # Initialize histograms

    hists = []

    hNTausPerEvent = TH1F('n_taus_per_event_' + decay_mode, 'Number of Taus Reconstructed Per Event (' + decay_mode + ')', 10, 0, 10)
    hists.append(hNTausPerEvent)
    hTypeCharged = TH1F('typeCharged_' + decay_mode, 'Type of Charged PFOs (' + decay_mode + ')', 5, 10, 15)
    hists.append(hTypeCharged)
    hNCharged = TH1F('nCharged_' + decay_mode, 'Number of Charged PFOs (' + decay_mode + ')', 10, 0, 10)
    hists.append(hNCharged)
    hD0Charged_11 = TH1F('d0Charged_pdg11_' + decay_mode, 'Impact Parameter of Charged PFOs (PDG 11, ' + decay_mode + ')', 100, 0, 1)
    hists.append(hD0Charged_11)
    hD0Charged_13 = TH1F('d0Charged_pdg13_' + decay_mode, 'Impact Parameter of Charged PFOs (PDG 13, ' + decay_mode + ')', 100, 0, 1)
    hists.append(hD0Charged_13)
    hD0Charged_211 = TH1F('d0Charged_pdg211_' + decay_mode, 'Impact Parameter of Charged PFOs (PDG 211, ' + decay_mode + ')', 100, 0, 1)
    hists.append(hD0Charged_211)
    hNTracks = TH1F('ntrks_' + decay_mode, 'Number of Reconstructed Tracks (' + decay_mode + ')', 7, 0, 7)
    hNTracks.GetXaxis().SetTitle('Number of Reconstructed Tracks')
    hists.append(hNTracks)
    hNPFOs = TH1F('nPFOs_' + decay_mode, 'Number of Reconstructed PFOs (' + decay_mode + ')', 35, 0, 35)
    hists.append(hNPFOs)
    hTypePFOs = TH1F('typePFOs_' + decay_mode, 'Type of Reconstructed PFOs (' + decay_mode + ')', 1100, 0, 2200)
    hists.append(hTypePFOs)
    hEnergyPFOs_11 = TH1F('energyPFOs_pdg11_' + decay_mode, 'Energy of Reconstructed PFOs (PDG 11, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hEnergyPFOs_11)
    hPtPFOs_11 = TH1F('ptPFOs_pdg11_' + decay_mode, 'Transverse Momentum of Reconstructed PFOs (PDG 11, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hPtPFOs_11)
    hEnergyRatioPFOs_11 = TH1F('energyRatioPFOs_pdg11_' + decay_mode, 'Energy Ratio of Reconstructed PFOs (PDG 11, ' + decay_mode + ')', 50, 0, 1)
    hists.append(hEnergyRatioPFOs_11)
    # hNHitsPFOs_11 = TH1F('nHitsPFOs_pdg11_' + decay_mode, 'Number of Detector Hits of Reconstructed PFOs (PDG 11, '+ decay_mode + ')', 250, 0, 500)
    # hists.append(hNHitsPFOs_11)
    hEnergyPFOs_13 = TH1F('energyPFOs_pdg13_' + decay_mode, 'Energy of Reconstructed PFOs (PDG 13, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hEnergyPFOs_13)
    hPtPFOs_13 = TH1F('ptPFOs_pdg13_' + decay_mode, 'Transverse Momentum of Reconstructed PFOs (PDG 13, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hPtPFOs_13)
    hEnergyRatioPFOs_13 = TH1F('energyRatioPFOs_pdg13_' + decay_mode, 'Energy Ratio of Reconstructed PFOs (PDG 13, ' + decay_mode + ')', 50, 0, 1)
    hists.append(hEnergyRatioPFOs_13)
    # hNHitsPFOs_13 = TH1F('nHitsPFOs_pdg13_' + decay_mode, 'Number of Detector Hits of Reconstructed PFOs (PDG 13, '+ decay_mode + ')', 250, 0, 500)
    # hists.append(hNHitsPFOs_13)
    hEnergyPFOs_22 = TH1F('energyPFOs_pdg22_' + decay_mode, 'Energy of Reconstructed PFOs (PDG 22, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hEnergyPFOs_22)
    hPtPFOs_22 = TH1F('ptPFOs_pdg22_' + decay_mode, 'Transverse Momentum of Reconstructed PFOs (PDG 22, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hPtPFOs_22)
    hEnergyRatioPFOs_22 = TH1F('energyRatioPFOs_pdg22_' + decay_mode, 'Energy Ratio of Reconstructed PFOs (PDG 22, ' + decay_mode + ')', 50, 0, 1)
    hists.append(hEnergyRatioPFOs_22)
    # hNHitsPFOs_22 = TH1F('nHitsPFOs_pdg22_' + decay_mode, 'Number of Detector Hits of Reconstructed PFOs (PDG 22, '+ decay_mode + ')', 250, 0, 500)
    # hists.append(hNHitsPFOs_22)
    hEnergyPFOs_211 = TH1F('energyPFOs_pdg211_' + decay_mode, 'Energy of Reconstructed PFOs (PDG 211, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hEnergyPFOs_211)
    hPtPFOs_211 = TH1F('ptPFOs_pdg211_' + decay_mode, 'Transverse Momentum of Reconstructed PFOs (PDG 211, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hPtPFOs_211)
    hEnergyRatioPFOs_211 = TH1F('energyRatioPFOs_pdg211_' + decay_mode, 'Energy Ratio of Reconstructed PFOs (PDG 211, ' + decay_mode + ')', 50, 0, 1)
    hists.append(hEnergyRatioPFOs_211)
    # hNHitsPFOs_211 = TH1F('nHitsPFOs_pdg211_' + decay_mode, 'Number of Detector Hits of Reconstructed PFOs (PDG 211, '+ decay_mode + ')', 250, 0, 500)
    # hists.append(hNHitsPFOs_211)
    hEnergyPFOs_2112 = TH1F('energyPFOs_pdg2112_' + decay_mode, 'Energy of Reconstructed PFOs (PDG 2112, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hEnergyPFOs_2112)
    hPtPFOs_2112 = TH1F('ptPFOs_pdg2112_' + decay_mode, 'Transverse Momentum of Reconstructed PFOs (PDG 2112, ' + decay_mode + ')', 50, 0, 10)
    hists.append(hPtPFOs_2112)
    hEnergyRatioPFOs_2112 = TH1F('energyRatioPFOs_pdg2112_' + decay_mode, 'Energy Ratio of Reconstructed PFOs (PDG 2112, ' + decay_mode + ')', 50, 0, 1)
    hists.append(hEnergyRatioPFOs_2112)
    # hNHitsPFOs_2112 = TH1F('nHitsPFOs_pdg2112_' + decay_mode, 'Number of Detector Hits of Reconstructed PFOs (PDG 2112, '+ decay_mode + ')', 250, 0, 500)
    # hists.append(hNHitsPFOs_2112)
    
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

    # Keep track of reconstructed PFO types
    pfo_types = {}

    # Keep track of charged PFO types
    charged_types = {}

    # Keep track of number of tracks
    n_trks = {}

    # Store visible mc tau properties
    E_vis = []
    px_vis = []
    py_vis = []
    pz_vis = []

    # Store number of reconstructed PFOs
    n_reco_PFOs = []

    # Store total number of taus that passed
    n_taus_passed = 0
    
    # Store total number of taus that failed number of charged tracks selection
    n_taus_nchargedtrks = 0
    
    # Store total number of taus that failed merge
    n_taus_merge = 0
    
    # Store total number of taus that failed number of particles selection
    n_taus_nparticles = 0

    types = []
    
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
            # tauRecoLink = event.getCollection('TauRecLink_PFO')
            # recoMCLink = event.getCollection('RecoMCTruthLink')

            # Instantiate relation navigators to parse tauReco and RecoMC links
            # relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
            # relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

            # Store number of reconstructed taus per event
            n_taus_per_event = 0

            for pfo in pfos:
              pfo_type = abs(pfo.getType())
              if pfo_type not in types:
                  types.append(pfo_type)
            
            # Loop through mcParticles
            for mcParticle in mcParticles:

                # Tag mcTaus
                if (abs(mcParticle.getPDG()) != 15):
                    continue
                    
                # Get true decay mode
                decay_mode_true = getDecayMode(mcParticle)

                # Tag desired decay mode
                if decay_mode_true == decay_mode_nums[n]:

                    # Store visible properties
                    visible_properties = getVisibleProperties(mcParticle)
                    E_vis.append(visible_properties[0])
                    px_vis.append(visible_properties[1])
                    py_vis.append(visible_properties[2])
                    pz_vis.append(visible_properties[3])
                        
                    reco_pfos = 0
                    
                    # Loop through rec taus which passed
                    for tau in taus:
                        n_taus_passed += 1
                        n_taus_per_event += 1

                        # Get reco tau energy
                        tau_energy = tau.getEnergy()

                        pfos = tau.getParticles()
                        npfos = len(pfos)
                        hNPFOs.Fill(npfos)
                        reco_pfos += npfos
                        ntrks = 0
                        for pfo in pfos:
                            pfo_type = abs(pfo.getType())
                            E = pfo.getEnergy()
                            E_ratio = E/tau_energy
                            px = pfo.getMomentum()[0]
                            py = pfo.getMomentum()[1]
                            pz = pfo.getMomentum()[2]
                            pt = math.sqrt(px**2 + py**2)
                            p = math.sqrt(pt**2 + pz**2)

                            if str(pfo_type) in pfo_types:
                                pfo_types[str(pfo_type)] += 1
                            else:
                                pfo_types[str(pfo_type)] = 1
                            
                            if pfo_type == 11:
                                hEnergyPFOs_11.Fill(E)
                                hPtPFOs_11.Fill(pt)
                                hEnergyRatioPFOs_11.Fill(E_ratio)
                                # hNHitsPFOs_11.Fill(n_hits)

                            elif pfo_type == 13:
                                hEnergyPFOs_13.Fill(E)
                                hPtPFOs_13.Fill(pt)
                                hEnergyRatioPFOs_13.Fill(E_ratio)
                                # hNHitsPFOs_13.Fill(n_hits)
                                    
                            elif pfo_type == 22:
                                hEnergyPFOs_22.Fill(E)
                                hPtPFOs_22.Fill(pt)
                                hEnergyRatioPFOs_22.Fill(E_ratio)
                                # hNHitsPFOs_22.Fill(n_hits)
                                    
                            elif pfo_type == 211:
                                hEnergyPFOs_211.Fill(E)
                                hPtPFOs_211.Fill(pt)
                                hEnergyRatioPFOs_211.Fill(E_ratio)
                                # hNHitsPFOs_211.Fill(n_hits)
                                    
                            elif pfo_type == 2112:
                                hEnergyPFOs_2112.Fill(E)
                                hPtPFOs_2112.Fill(pt)
                                hEnergyRatioPFOs_2112.Fill(E_ratio)
                                # hNHitsPFOs_2112.Fill(n_hits)

                            charge = pfo.getCharge()
                            if charge != 0:
                                tracks = pfo.getTracks()
                                ntrks += len(tracks)
                                for track in tracks:
                                    d0 = track.getD0()

                                    if pfo_type == 11:
                                        hD0Charged_11.Fill(d0)
                                    elif pfo_type == 13:
                                        hD0Charged_13.Fill(d0)
                                    elif pfo_type == 211:
                                        hD0Charged_211.Fill(d0)

                        hNTracks.Fill(ntrks)

                        if str(ntrks) in n_trks:
                            n_trks[str(ntrks)] += 1
                        else:
                            n_trks[str(ntrks)] = 1
                                
                    # Loop through reco taus which failed number of charged tracks
                    for tau_nchargedtrks in taus_nchargedtrks:

                        n_taus_nchargedtrks += 1
                        n_taus_per_event += 1

                        # Get reco tau energy
                        tau_energy = tau_nchargedtrks.getEnergy()
                        
                        pfos = tau_nchargedtrks.getParticles()
                        
                        ncharged = 0
                        for pfo in pfos:
                            if pfo.getCharge() == 0:
                                continue
                            ncharged += 1

                            charged_type = abs(pfo.getType())
                            hTypeCharged.Fill(charged_type)

                            if str(charged_type) in charged_types:
                                charged_types[str(charged_type)] += 1
                            else:
                                charged_types[str(charged_type)] = 1

                        hNCharged.Fill(ncharged)

                        npfos = len(pfos)
                        hNPFOs.Fill(npfos)
                        reco_pfos += npfos
                        ntrks = 0
                        for pfo in pfos:
                            pfo_type = abs(pfo.getType())
                            E = pfo.getEnergy()
                            E_ratio = E/tau_energy
                            px = pfo.getMomentum()[0]
                            py = pfo.getMomentum()[1]
                            pz = pfo.getMomentum()[2]
                            pt = math.sqrt(px**2 + py**2)
                            p = math.sqrt(pt**2 + pz**2)

                            if str(pfo_type) in pfo_types:
                                pfo_types[str(pfo_type)] += 1
                            else:
                                pfo_types[str(pfo_type)] = 1

                            if pfo_type == 11:
                                hEnergyPFOs_11.Fill(E)
                                hPtPFOs_11.Fill(pt)
                                hEnergyRatioPFOs_11.Fill(E_ratio)
                                # hNHitsPFOs_11.Fill(n_hit)

                            elif pfo_type == 13:
                                hEnergyPFOs_13.Fill(E)
                                hPtPFOs_13.Fill(pt)
                                hEnergyRatioPFOs_13.Fill(E_ratio)
                                # hNHitsPFOs_13.Fill(n_hit)
                                
                            elif pfo_type == 22:
                                hEnergyPFOs_22.Fill(E)
                                hPtPFOs_22.Fill(pt)
                                hEnergyRatioPFOs_22.Fill(E_ratio)
                                # hNHitsPFOs_22.Fill(n_hit)
                                
                            elif pfo_type == 211:
                                hEnergyPFOs_211.Fill(E)
                                hPtPFOs_211.Fill(pt)
                                hEnergyRatioPFOs_211.Fill(E_ratio)
                                # hNHitsPFOs_211.Fill(n_hit)
                                
                            elif pfo_type == 2112:
                                hEnergyPFOs_2112.Fill(E)
                                hPtPFOs_2112.Fill(pt)
                                hEnergyRatioPFOs_2112.Fill(E_ratio)
                                # hNHitsPFOs_2112.Fill(n_hit)
                                
                            charge = pfo.getCharge()
                            if charge != 0:
                                tracks = pfo.getTracks()
                                ntrks += len(tracks)
                                for track in tracks:
                                    d0 = track.getD0()

                                    if pfo_type == 11:
                                        hD0Charged_11.Fill(d0)
                                    elif pfo_type == 13:
                                        hD0Charged_13.Fill(d0)
                                    elif pfo_type == 211:
                                        hD0Charged_211.Fill(d0)

                        hNTracks.Fill(ntrks)

                        if str(ntrks) in n_trks:
                            n_trks[str(ntrks)] += 1
                        else:
                            n_trks[str(ntrks)] = 1
                        
                    # Loop through reco taus which failed merge
                    for tau_merge in taus_merge:

                        n_taus_merge += 1
                        n_taus_per_event += 1
                        
                        # Get reco tau energy
                        tau_energy = tau_merge.getEnergy()
                        
                        pfos = tau_merge.getParticles()
                        npfos = len(pfos)
                        hNPFOs.Fill(npfos)
                        reco_pfos += npfos
                        ntrks = 0
                        for pfo in pfos:
                            pfo_type = abs(pfo.getType())
                            E = pfo.getEnergy()
                            E_ratio = E/tau_energy
                            px = pfo.getMomentum()[0]
                            py = pfo.getMomentum()[1]
                            pz = pfo.getMomentum()[2]
                            pt = math.sqrt(px**2 + py**2)
                            p = math.sqrt(pt**2 + pz**2)
                            
                            if str(pfo_type) in pfo_types:
                                pfo_types[str(pfo_type)] += 1
                            else:
                                pfo_types[str(pfo_type)] = 1
                            
                            if pfo_type == 11:
                                hEnergyPFOs_11.Fill(E)
                                hPtPFOs_11.Fill(pt)
                                hEnergyRatioPFOs_11.Fill(E_ratio)
                                # hNHitsPFOs_11.Fill(n_hits)

                            elif pfo_type == 13:
                                hEnergyPFOs_13.Fill(E)
                                hPtPFOs_13.Fill(pt)
                                hEnergyRatioPFOs_13.Fill(E_ratio)
                                # hNHitsPFOs_13.Fill(n_hits)
                                
                            elif pfo_type == 22:
                                hEnergyPFOs_22.Fill(E)
                                hPtPFOs_22.Fill(pt)
                                hEnergyRatioPFOs_22.Fill(E_ratio)
                                # hNHitsPFOs_22.Fill(n_hits)
                                
                            elif pfo_type == 211:
                                hEnergyPFOs_211.Fill(E)
                                hPtPFOs_211.Fill(pt)
                                hEnergyRatioPFOs_211.Fill(E_ratio)
                                # hNHitsPFOs_211.Fill(n_hits)
                                
                            elif pfo_type == 2112:
                                hEnergyPFOs_2112.Fill(E)
                                hPtPFOs_2112.Fill(pt)
                                hEnergyRatioPFOs_2112.Fill(E_ratio)
                                # hNHitsPFOs_2112.Fill(n_hits)
                                
                            charge = pfo.getCharge()
                            if charge != 0:
                                tracks = pfo.getTracks()
                                ntrks += len(tracks)
                                for track in tracks:
                                    d0 = track.getD0()

                                    if pfo_type == 11:
                                        hD0Charged_11.Fill(d0)
                                    elif pfo_type == 13:
                                        hD0Charged_13.Fill(d0)
                                    elif pfo_type == 211:
                                        hD0Charged_211.Fill(d0)

                        hNTracks.Fill(ntrks)

                        if str(ntrks) in n_trks:
                            n_trks[str(ntrks)] += 1
                        else:
                            n_trks[str(ntrks)] = 1
                                
                    # Loop through reco taus which failed number of particles
                    for tau_nparticles in taus_nparticles:

                        n_taus_nparticles += 1
                        n_taus_per_event += 1
                        
                        # Get reco tau energy
                        tau_energy = tau_nparticles.getEnergy()
                        
                        pfos = tau_nparticles.getParticles()

                        npfos = len(pfos)
                        hNPFOs.Fill(npfos)

                        reco_pfos += npfos
                        ntrks = 0
                        for pfo in pfos:
                            pfo_type = abs(pfo.getType())
                            hTypePFOs.Fill(pfo_type)

                            if str(pfo_type) in pfo_types:
                                pfo_types[str(pfo_type)] += 1
                            else:
                                pfo_types[str(pfo_type)] = 1

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
                                # hNHitsPFOs_11.Fill(n_hits)

                            elif pfo_type == 13:
                                hEnergyPFOs_13.Fill(E)
                                hPtPFOs_13.Fill(pt)
                                hEnergyRatioPFOs_13.Fill(E_ratio)
                                # hNHitsPFOs_13.Fill(n_hits)
                                
                            elif pfo_type == 22:
                                hEnergyPFOs_22.Fill(E)
                                hPtPFOs_22.Fill(pt)
                                hEnergyRatioPFOs_22.Fill(E_ratio)
                                # hNHitsPFOs_22.Fill(n_hits)
                                
                            elif pfo_type == 211:
                                hEnergyPFOs_211.Fill(E)
                                hPtPFOs_211.Fill(pt)
                                hEnergyRatioPFOs_211.Fill(E_ratio)
                                # hNHitsPFOs_211.Fill(n_hits)
                                
                            elif pfo_type == 2112:
                                hEnergyPFOs_2112.Fill(E)
                                hPtPFOs_2112.Fill(pt)
                                hEnergyRatioPFOs_2112.Fill(E_ratio)
                                # hNHitsPFOs_2112.Fill(n_hits)
                                
                            charge = pfo.getCharge()
                            if charge != 0:
                                tracks = pfo.getTracks()
                                ntrks += len(tracks)
                                for track in tracks:
                                    d0 = track.getD0()

                                    if pfo_type == 11:
                                        hD0Charged_11.Fill(d0)
                                    elif pfo_type == 13:
                                        hD0Charged_13.Fill(d0)
                                    elif pfo_type == 211:
                                        hD0Charged_211.Fill(d0)

                        hNTracks.Fill(ntrks)

                        if str(ntrks) in n_trks:
                            n_trks[str(ntrks)] += 1
                        else:
                            n_trks[str(ntrks)] = 1
                                        
                    n_reco_PFOs.append(reco_pfos)

            hNTausPerEvent.Fill(n_taus_per_event)
                    
        # Close file
        reader.close()

    for pfo_type in types:
        print(f'PFO Type: {pfo_type}')
        
    # Calculate remaining visible properties
    pt_vis = np.sqrt(np.power(px_vis, 2) + np.power(py_vis, 2))
    p_vis = np.sqrt(np.power(pt_vis, 2) + np.power(pz_vis, 2))
    theta_vis = np.arccos(pz_vis/p_vis)
    phi_vis = np.arccos(px_vis/pt_vis)

    # Plot nPFOs versus visible observables
    plt.figure()

    plt.clf()
    plt.scatter(E_vis, n_reco_PFOs, color='blue', s=8)
    plt.axhline(y=10, linestyle='--', color='red')
    plt.xlabel(r'True $\tau^-$ Visible Energy [GeV]')
    plt.ylabel('Number of PFOs Reconstructed')
    plt.title(decay_mode + ' Number of PFOs vs Visible E')
    plt.savefig(decay_mode + '_nPFOs_E_vis.png')

    plt.clf()
    plt.scatter(np.absolute(px_vis), n_reco_PFOs, color='blue', s=8)
    plt.axhline(y=10, linestyle='--', color='red')
    plt.xlabel(r'True $\tau^-$ Visible $p_x$ [GeV/c]')
    plt.ylabel('Number of PFOs Reconstructed')
    plt.title(decay_mode + r' Number of PFOs vs Visible $p_x$')
    plt.savefig(decay_mode + '_nPFOs_px_vis.png')

    plt.clf()
    plt.scatter(np.absolute(py_vis), n_reco_PFOs, color='blue', s=8)
    plt.axhline(y=10, linestyle='--', color='red')
    plt.xlabel(r'True $\tau^-$ Visible $p_y$ [GeV/c]')
    plt.ylabel('Number of PFOs Reconstructed')
    plt.title(decay_mode + r' Number of PFOs vs Visible $p_y$')
    plt.savefig(decay_mode + '_nPFOs_py_vis.png')

    plt.clf()
    plt.scatter(np.absolute(pz_vis), n_reco_PFOs, color='blue', s=8)
    plt.axhline(y=10, linestyle='--', color='red')
    plt.xlabel(r'True $\tau^-$ Visible $p_z$ [GeV/c]')
    plt.ylabel('Number of PFOs Reconstructed')
    plt.title(decay_mode + r' Number of PFOs vs Visible $p_z$')
    plt.savefig(decay_mode + '_nPFOs_pz_vis.png')

    plt.clf()
    plt.scatter(np.absolute(pt_vis), n_reco_PFOs, color='blue', s=8)
    plt.axhline(y=10, linestyle='--', color='red')
    plt.xlabel(r'True $\tau^-$ Visible $p_t$ [GeV/c]')
    plt.ylabel('Number of PFOs Reconstructed')
    plt.title(decay_mode + r' Number of PFOs vs Visible $p_t$')
    plt.savefig(decay_mode + '_nPFOs_pt_vis.png')

    plt.clf()
    plt.scatter(np.absolute(p_vis), n_reco_PFOs, color='blue', s=8)
    plt.axhline(y=10, linestyle='--', color='red')
    plt.xlabel(r'True $\tau^-$ Visible $p$ [GeV/c]')
    plt.ylabel('Number of PFOs Reconstructed')
    plt.title(decay_mode + r' Number of PFOs vs Visible $p$')
    plt.savefig(decay_mode + '_nPFOs_p_vis.png')

    plt.clf()
    plt.scatter(np.absolute(theta_vis), n_reco_PFOs, color='blue', s=8)
    plt.axhline(y=10, linestyle='--', color='red')
    plt.xlabel(r'True $\tau^-$ Visible $\theta$ [rad]')
    plt.ylabel('Number of PFOs Reconstructed')
    plt.title(decay_mode + r' Number of PFOs vs Visible $\theta$')
    plt.savefig(decay_mode + '_nPFOs_theta_vis.png')

    plt.clf()
    plt.scatter(np.absolute(phi_vis), n_reco_PFOs, color='blue', s=8)
    plt.axhline(y=10, linestyle='--', color='red')
    plt.xlabel(r'True $\tau^-$ Visible $\phi$ [rad]')
    plt.ylabel('Number of PFOs Reconstructed')
    plt.title(decay_mode + r' Number of PFOs vs Visible $\phi$')
    plt.savefig(decay_mode + '_nPFOs_phi_vis.png')
    
    # Write to output file
    output_file = TFile(decay_mode + '_ana.root', 'RECREATE')
    for hist in hists:
        hist.Write()
    output_file.Close()

    # Draw hists and save as PNG
    for hist in hists:
        filename = hist.GetName() + '.png'
        canvas = TCanvas()
        hist.Draw()
        canvas.SaveAs(filename)

    with open(decay_mode + '_pfo_types.txt', 'w') as file:
        for key, value in pfo_types.items():
            print(f'Type: {key}, Total Number: {value}', file=file)

    with open(decay_mode + '_charged_types.txt', 'w') as file:
        for key, value in charged_types.items():
            print(f'Type: {key}, Total Number: {value}', file=file)

    with open(decay_mode + '_cutflow.txt', 'w') as file:
        print(f'Number of taus that passed: {n_taus_passed}', file=file)
        print(f'Number of taus that failed number of charged tracks: {n_taus_nchargedtrks}', file=file)
        print(f'Number of taus that failed merge: {n_taus_merge}', file=file)
        print(f'Number of taus that failed number of particles: {n_taus_nparticles}', file=file)

    with open(decay_mode + '_ntrks.txt', 'w') as file:
        for key, value in n_trks.items():
            print(f'Number of Tracks: {key}, Total Numebr: {value}', file=file)
