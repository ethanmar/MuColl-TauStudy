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

    hNPFOs = TH1F('nPFOs_' + decay_mode, 'Number of Reconstructed PFOs (' + decay_mode + ')', 50, 0, 50)
    hists.append(hNPFOs)
    hTypePFOs = TH1F('typePFOs_' + decay_mode, 'Type of Reconstructed PFOs (' + decay_mode + ')', 1100, 0, 2200)
    hists.append(hTypePFOs)

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

    # Store visible mc tau properties
    E_vis = []
    px_vis = []
    py_vis = []
    pz_vis = []

    # Store number of reconstructed PFOs
    n_reco_PFOs = []

    # Open input file(s)
    for file in to_process:
        reader = IOIMPL.LCFactory.getInstance().createLCReader()
        reader.open(file)

        # Loop through events
        for ievt, event in enumerate(reader):

            # Get collections
            # taus = event.getCollection('TauRec_PFO')
            pfos = event.getCollection('PandoraPFOs')
            mcParticles = event.getCollection('MCParticle')
            # tauRecoLink = event.getCollection('TauRecLink_PFO')
            # recoMCLink = event.getCollection('RecoMCTruthLink')

            # Instantiate relation navigators to parse tauReco and RecoMC links
            # relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
            # relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

            # Loop through mcParticles
            for mcParticle in mcParticles:

                # Tag mcTaus
                if (abs(mcParticle.getPDG()) == 15):
                    
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
                        
                        # Fill number of reconstructed PFOs
                        nPFOs = len(pfos)
                        n_reco_PFOs.append(nPFOs)
                        hNPFOs.Fill(nPFOs)

                        # if nPFOs >= 10:
                            # print(f'{decay_mode} NUMBER OF PFOS EXCEEDS OR EQUALS 10!!! Number of PFOs in event: {nPFOs}. Run number: {event.getRunNumber()}. Event number: {event.getEventNumber()}')
                        
                        # Loop over reconstructed PFOs
                        for pfo in pfos:

                            # Fill types of reconstructed PFOs
                            pfo_type = abs(pfo.getType())
                            hTypePFOs.Fill(pfo_type)

                            if str(pfo_type) in pfo_types:
                                pfo_types[str(pfo_type)] += 1
                            else:
                                pfo_types[str(pfo_type)] = 1

        # Close file
        reader.close()

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
