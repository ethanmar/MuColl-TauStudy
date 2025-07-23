from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

from tau_mc_link import getLinkedMCTau, getDecayMode

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
        if (pdg == 12 or pdg == 14 or pdg == 16):
            continue

        # Sum visible observables
        E_vis += mcDaughter.getEnergy()
        px_vis += mcDaughter.getMomentum()[0]
        py_vis += mcDaughter.getMomentum()[1]
        pz_vis += mcDaughter.getMomentum()[2]

    vis_properties = [E_vis, px_vis, py_vis, pz_vis]
    return vis_properties

def getPt(vis_props):
    return math.sqrt(vis_props[1]**2 + vis_props[2]**2)

def getP(vis_props):
    return math.sqrt(vis_props[1]**2 + vis_props[2]**2 + vis_props[3]**2)

def getTheta(vis_props):
    p = getP(vis_props)
    pz = vis_props[3]
    return math.acos(pz/p)

def getPhi(vis_props):
    px = vis_props[1]
    pt = getPt(vis_props)
    return math.acos(px/pt)

def getDetectorRegion(theta):
    if (1.0 < theta and theta < 2.0):
        return "Central Barrel"
    elif ((0.577 < theta and theta < 1.0) or (2.0 < theta and theta < 2.56)):
        return "Transition"
    elif (theta < 0.577 or theta > 2.56):
        return "Endcap"
    else:
        return "Other"

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='pi_eff.root')

args = parser.parse_args()

# Initialize histograms
hists = []

# 1P0N Pions
hPiPtTrue1P0N = TH1F('1p0n_pi_true_pT', 'True 1P0N Pion Pt', 10, 0, 320)
hists.append(hPiPtTrue1P0N)

hPiPhiTrue1P0N = TH1F('1p0n_pi_true_phi', 'True 1P0N Pion Phi', 10, 0, math.pi)
hists.append(hPiPhiTrue1P0N)

hPiThetaTrue1P0N = TH1F('1p0n_pi_true_theta', 'True 1P0N Pion Theta', 10, 0, math.pi)
hists.append(hPiThetaTrue1P0N)

hPiPtMatched1P0N = TH1F('1p0n_pi_matched_pT', 'Matched 1P0N Pion Pt', 10, 0, 320)
hists.append(hPiPtMatched1P0N)

hPiPhiMatched1P0N = TH1F('1p0n_pi_matched_phi', 'Matched 1P0N Pion Phi', 10, 0, math.pi)
hists.append(hPiPhiMatched1P0N)

hPiThetaMatched1P0N = TH1F('1p0n_pi_matched_theta', 'Matched 1P0N Pion Theta', 10, 0, math.pi)
hists.append(hPiThetaMatched1P0N)

hPiPtTrue1P0N_Barrel = TH1F('1p0n_pi_true_pT_barrel', 'True 1P0N Pion Pt (Barrel)', 10, 0, 320)
hists.append(hPiPtTrue1P0N_Barrel)

hPiPhiTrue1P0N_Barrel = TH1F('1p0n_pi_true_phi_barrel', 'True 1P0N Pion Phi (Barrel)', 10, 0, math.pi)
hists.append(hPiPhiTrue1P0N_Barrel)

hPiPtMatched1P0N_Barrel = TH1F('1p0n_pi_matched_pT_barrel', 'Matched 1P0N Pion Pt (Barrel)', 10, 0, 320)
hists.append(hPiPtMatched1P0N_Barrel)

hPiPhiMatched1P0N_Barrel = TH1F('1p0n_pi_matched_phi_barrel', 'Matched 1P0N Pion Phi (Barrel)', 10, 0, math.pi)
hists.append(hPiPhiMatched1P0N_Barrel)

hPiPtTrue1P0N_CentBarrel = TH1F('1p0n_pi_true_pT_centbarrel', 'True 1P0N Pion Pt (Central Barrel)', 10, 0, 320)
hists.append(hPiPtTrue1P0N_CentBarrel)

hPiPhiTrue1P0N_CentBarrel = TH1F('1p0n_pi_true_phi_centbarrel', 'True 1P0N Pion Phi (Central Barrel)', 10, 0, math.pi)
hists.append(hPiPhiTrue1P0N_CentBarrel)

hPiPtMatched1P0N_CentBarrel = TH1F('1p0n_pi_matched_pT_centbarrel', 'Matched 1P0N Pion Pt (Central Barrel)', 10, 0, 320)
hists.append(hPiPtMatched1P0N_CentBarrel)

hPiPhiMatched1P0N_CentBarrel = TH1F('1p0n_pi_matched_phi_centbarrel', 'Matched 1P0N Pion Phi (Central Barrel)', 10, 0, math.pi)
hists.append(hPiPhiMatched1P0N_CentBarrel)

hPiPtTrue1P0N_Endcap = TH1F('1p0n_pi_true_pT_endcap', 'True 1P0N Pion Pt (Endcap)', 10, 0, 320)
hists.append(hPiPtTrue1P0N_Endcap)

hPiPhiTrue1P0N_Endcap = TH1F('1p0n_pi_true_phi_endcap', 'True 1P0N Pion Phi (Endcap)', 10, 0, math.pi)
hists.append(hPiPhiTrue1P0N_Endcap)

hPiPtMatched1P0N_Endcap = TH1F('1p0n_pi_matched_pT_endcap', 'Matched 1P0N Pion Pt (Endcap)', 10, 0, 320)
hists.append(hPiPtMatched1P0N_Endcap)

hPiPhiMatched1P0N_Endcap = TH1F('1p0n_pi_matched_phi_endcap', 'Matched 1P0N Pion Phi (Endcap)', 10, 0, math.pi)
hists.append(hPiPhiMatched1P0N_Endcap)

hPiPtTrue1P0N_Transition = TH1F('1p0n_pi_true_pT_transition', 'True 1P0N Pion Pt (Transition)', 10, 0, 320)
hists.append(hPiPtTrue1P0N_Transition)

hPiPhiTrue1P0N_Transition = TH1F('1p0n_pi_true_phi_transition', 'True 1P0N Pion Phi (Transition)', 10, 0, math.pi)
hists.append(hPiPhiTrue1P0N_Transition)

hPiPtMatched1P0N_Transition = TH1F('1p0n_pi_matched_pT_transition', 'Matched 1P0N Pion Pt (Transition)', 10, 0, 320)
hists.append(hPiPtMatched1P0N_Transition)

hPiPhiMatched1P0N_Transition = TH1F('1p0n_pi_matched_phi_transition', 'Matched 1P0N Pion Phi (Transition)', 10, 0, math.pi)
hists.append(hPiPhiMatched1P0N_Transition)

# 3P0N Pions
hPiPtTrue3P0N = TH1F('3p0n_pi_true_pT', 'True 3P0N Pion Pt', 10, 0, 320)
hists.append(hPiPtTrue3P0N)

hPiPhiTrue3P0N = TH1F('3p0n_pi_true_phi', 'True 3P0N Pion Phi', 10, 0, math.pi)
hists.append(hPiPhiTrue3P0N)

hPiThetaTrue3P0N = TH1F('3p0n_pi_true_theta', 'True 3P0N Pion Theta', 10, 0, math.pi)
hists.append(hPiThetaTrue3P0N)

hPiPtMatched3P0N = TH1F('3p0n_pi_matched_pT', 'Matched 3P0N Pion Pt', 10, 0, 320)
hists.append(hPiPtMatched3P0N)

hPiPhiMatched3P0N = TH1F('3p0n_pi_matched_phi', 'Matched 3P0N Pion Phi', 10, 0, math.pi)
hists.append(hPiPhiMatched3P0N)

hPiThetaMatched3P0N = TH1F('3p0n_pi_matched_theta', 'Matched 3P0N Pion Theta', 10, 0, math.pi)
hists.append(hPiThetaMatched3P0N)

hPiPtTrue3P0N_Barrel = TH1F('3p0n_pi_true_pT_barrel', 'True 3P0N Pion Pt (Barrel)', 10, 0, 320)
hists.append(hPiPtTrue3P0N_Barrel)

hPiPhiTrue3P0N_Barrel = TH1F('3p0n_pi_true_phi_barrel', 'True 3P0N Pion Phi (Barrel)', 10, 0, math.pi)
hists.append(hPiPhiTrue3P0N_Barrel)

hPiPtMatched3P0N_Barrel = TH1F('3p0n_pi_matched_pT_barrel', 'Matched 3P0N Pion Pt (Barrel)', 10, 0, 320)
hists.append(hPiPtMatched3P0N_Barrel)

hPiPhiMatched3P0N_Barrel = TH1F('3p0n_pi_matched_phi_barrel', 'Matched 3P0N Pion Phi (Barrel)', 10, 0, math.pi)
hists.append(hPiPhiMatched3P0N_Barrel)

hPiPtTrue3P0N_CentBarrel = TH1F('3p0n_pi_true_pT_centbarrel', 'True 3P0N Pion Pt (Central Barrel)', 10, 0, 320)
hists.append(hPiPtTrue3P0N_CentBarrel)

hPiPhiTrue3P0N_CentBarrel = TH1F('3p0n_pi_true_phi_centbarrel', 'True 3P0N Pion Phi (Central Barrel)', 10, 0, math.pi)
hists.append(hPiPhiTrue3P0N_CentBarrel)

hPiPtMatched3P0N_CentBarrel = TH1F('3p0n_pi_matched_pT_centbarrel', 'Matched 3P0N Pion Pt (Central Barrel)', 10, 0, 320)
hists.append(hPiPtMatched3P0N_CentBarrel)

hPiPhiMatched3P0N_CentBarrel = TH1F('3p0n_pi_matched_phi_centbarrel', 'Matched 3P0N Pion Phi (Central Barrel)', 10, 0, math.pi)
hists.append(hPiPhiMatched3P0N_CentBarrel)

hPiPtTrue3P0N_Endcap = TH1F('3p0n_pi_true_pT_endcap', 'True 3P0N Pion Pt (Endcap)', 10, 0, 320)
hists.append(hPiPtTrue3P0N_Endcap)

hPiPhiTrue3P0N_Endcap = TH1F('3p0n_pi_true_phi_endcap', 'True 3P0N Pion Phi (Endcap)', 10, 0, math.pi)
hists.append(hPiPhiTrue3P0N_Endcap)

hPiPtMatched3P0N_Endcap = TH1F('3p0n_pi_matched_pT_endcap', 'Matched 3P0N Pion Pt (Endcap)', 10, 0, 320)
hists.append(hPiPtMatched3P0N_Endcap)

hPiPhiMatched3P0N_Endcap = TH1F('3p0n_pi_matched_phi_endcap', 'Matched 3P0N Pion Phi (Endcap)', 10, 0, math.pi)
hists.append(hPiPhiMatched3P0N_Endcap)

hPiPtTrue3P0N_Transition = TH1F('3p0n_pi_true_pT_transition', 'True 3P0N Pion Pt (Transition)', 10, 0, 320)
hists.append(hPiPtTrue3P0N_Transition)

hPiPhiTrue3P0N_Transition = TH1F('3p0n_pi_true_phi_transition', 'True 3P0N Pion Phi (Transition)', 10, 0, math.pi)
hists.append(hPiPhiTrue3P0N_Transition)

hPiPtMatched3P0N_Transition = TH1F('3p0n_pi_matched_pT_transition', 'Matched 3P0N Pion Pt (Transition)', 10, 0, 320)
hists.append(hPiPtMatched3P0N_Transition)

hPiPhiMatched3P0N_Transition = TH1F('3p0n_pi_matched_phi_transition', 'Matched 3P0N Pion Phi (Transition)', 10, 0, math.pi)
hists.append(hPiPhiMatched3P0N_Transition)

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

n_1p0n_true = 0
n_1p0n_reco = 0
n_3p0n_true = 0
n_3p0n_reco = 0

n_pi_true = []
n_pi_reco = []

n_electrons = []
tau_pts = []

# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        taus = event.getCollection('RecoTaus')
        pfos = event.getCollection('PandoraPFOs')
        mc_particles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauPFOLink')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        # Instantiate relation navigators to parse tauReco and RecoMC links
        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        # Store visible tau properties for plotting
        vis_props = []
        vis_pt = 0
        vis_theta = 0
        vis_phi = 0

        # Store decay mode of event
        decay_mode = -1

        # True pi observables
        true_pi_pts = []
        true_pi_thetas = []
        true_pis = []
        true_pi_phis = []
        true_pi_charges = []
        true_pi_energies = []
        
        # Loop through mc_particles
        for mc_particle in mc_particles:

            # Tag taus
            if (abs(mc_particle.getPDG()) == 15):

                # Get decay mode
                decay_mode = getDecayMode(mc_particle)

                if (decay_mode == 0 or decay_mode == 4):

                    # Get visible properties
                    vis_props = getVisibleProperties(mc_particle)
                    vis_pt = getPt(vis_props)
                    vis_theta = getTheta(vis_props)
                    vis_phi = getPhi(vis_props)

                    if decay_mode == 4:
                        daughters = mc_particle.getDaughters()
                        for daughter in daughters:
                            if (abs(daughter.getPDG()) == 211):
                                mom = daughter.getMomentum()
                                pt = math.sqrt(mom[0]**2 + mom[1]**2)
                                p = math.sqrt(pt**2 + mom[2]**2)
                                pz = mom[2]
                                px = mom[0]
                                theta = math.acos(pz/p)
                                phi = math.acos(px/pt)
                                charge = daughter.getCharge()
                                energy = daughter.getEnergy()
                                true_pi_pts.append(pt)
                                true_pi_thetas.append(theta)
                                true_pis.append(daughter)
                                true_pi_phis.append(phi)
                                true_pi_charges.append(charge)
                                true_pi_energies.append(energy)

        if decay_mode == 0:
            n_pi_true.append(1)
            n_1p0n_true += 1
            hPiPtTrue1P0N.Fill(vis_pt)
            hPiThetaTrue1P0N.Fill(vis_theta)
            hPiPhiTrue1P0N.Fill(vis_phi)

            if (vis_theta > 0.70 and vis_theta < 2.45):
                hPiPtTrue1P0N_Barrel.Fill(vis_pt)
                hPiPhiTrue1P0N_Barrel.Fill(vis_phi)

            if (vis_theta > 1 and vis_theta < 2):
                hPiPtTrue1P0N_CentBarrel.Fill(vis_pt)
                hPiPhiTrue1P0N_CentBarrel.Fill(vis_phi)

            elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                hPiPtTrue1P0N_Transition.Fill(vis_pt)
                hPiPhiTrue1P0N_Transition.Fill(vis_phi)

            elif (vis_theta < 0.577 or vis_theta > 2.56):
                hPiPtTrue1P0N_Endcap.Fill(vis_pt)
                hPiPhiTrue1P0N_Endcap.Fill(vis_phi)

        elif decay_mode == 4:
            tau_pts.append(vis_pt)
            n_pi_true.append(3)
            n_3p0n_true += 1
            hPiPtTrue3P0N.Fill(vis_pt)
            hPiThetaTrue3P0N.Fill(vis_theta)
            hPiPhiTrue3P0N.Fill(vis_phi)

            if (vis_theta > 0.70 and vis_theta < 2.45):
                hPiPtTrue3P0N_Barrel.Fill(vis_pt)
                hPiPhiTrue3P0N_Barrel.Fill(vis_phi)

            if (vis_theta > 1 and vis_theta < 2):
                hPiPtTrue3P0N_CentBarrel.Fill(vis_pt)
                hPiPhiTrue3P0N_CentBarrel.Fill(vis_phi)

            elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                hPiPtTrue3P0N_Transition.Fill(vis_pt)
                hPiPhiTrue3P0N_Transition.Fill(vis_phi)

            elif (vis_theta < 0.577 or vis_theta > 2.56):
                hPiPtTrue3P0N_Endcap.Fill(vis_pt)
                hPiPhiTrue3P0N_Endcap.Fill(vis_phi)

        # Keep track of number of matched reco charged pion
        n_matched_pis_1p0n = 0
        n_matched_duplicate_pis_3p0n = 0
        n_matched_unique_pis_3p0n = 0
        n_matched_pis_all = 0
        
        # Keep track of true, matched 3p0n pis
        matched_true_3p0n_pis = []

        unique_pts = []
        unique_energies = []

        duplicate_pts = []
        duplicate_energies = []
        
        # Misidentified
        mis_id_types = []
        mis_id_pts = []
        mis_id_energies = []

        n_matched_electrons = 0
        
        # Loop through pfos
        for pfo in pfos:

            # Tag charged pions
            if (abs(pfo.getType()) == 211):

                # Get linked mc particles
                linked_mc_particles = relationNavigatorRecoMC.getRelatedToObjects(pfo)
                linked_weights = relationNavigatorRecoMC.getRelatedToWeights(pfo)

                # Get most weighted linked mc particle
                max_weight = 0
                max_weight_mc_particle = None
                for linked_weight, linked_mc_particle in zip(linked_weights, linked_mc_particles):
                    if linked_weight > max_weight:
                        max_weight = linked_weight
                        max_weight_mc_particle = linked_mc_particle

                # Ensure matched particle exists
                if max_weight_mc_particle is None:
                    continue
                        
                # Check if most weighted mc particle is a charged pion        
                if (abs(max_weight_mc_particle.getPDG()) == 211):

                    # Check if parent is a tau
                    parent_tau = None
                    parents = max_weight_mc_particle.getParents()
                    for parent in parents:
                        if (abs(parent.getPDG()) == 15):
                            parent_tau = parent

                    # Check decay mode of parent tau
                    if parent_tau is not None:

                        decay_mode = getDecayMode(parent_tau)

                        if decay_mode == 0:
                            n_matched_pis_1p0n += 1
                        elif decay_mode == 4:
                            if max_weight_mc_particle not in matched_true_3p0n_pis:
                                matched_true_3p0n_pis.append(max_weight_mc_particle)
                                n_matched_unique_pis_3p0n += 1
                                mom = pfo.getMomentum()
                                pt = math.sqrt(mom[0]**2 + mom[1]**2)
                                energy = pfo.getEnergy()
                                unique_pts.append(pt)
                                unique_energies.append(energy)
                            else:
                                n_matched_duplicate_pis_3p0n += 1
                                mom = pfo.getMomentum()
                                pt = math.sqrt(mom[0]**2 + mom[1]**2)
                                energy = pfo.getEnergy()
                                duplicate_pts.append(pt)
                                duplicate_energies.append(energy)
            else:

                # Get pfo type
                pfo_type = abs(pfo.getType())
                
                # Get linked mc particles
                linked_mc_particles = relationNavigatorRecoMC.getRelatedToObjects(pfo)
                linked_weights = relationNavigatorRecoMC.getRelatedToWeights(pfo)

                # Get most weighted linked mc particle
                max_weight = 0
                max_weight_mc_particle = None
                for linked_weight, linked_mc_particle in zip(linked_weights, linked_mc_particles):
                    if linked_weight > max_weight:
                        max_weight = linked_weight
                        max_weight_mc_particle = linked_mc_particle

                # Ensure matched particle exists
                if max_weight_mc_particle is None:
                    continue
                        
                # Check if most weighted mc particle is a charged pion        
                if (abs(max_weight_mc_particle.getPDG()) == 211):

                    # Check if parent is a tau
                    parent_tau = None
                    parents = max_weight_mc_particle.getParents()
                    for parent in parents:
                        if (abs(parent.getPDG()) == 15):
                            parent_tau = parent

                    # Check decay mode of parent tau
                    if parent_tau is not None:

                        decay_mode = getDecayMode(parent_tau)

                        if decay_mode == 4:
                            mis_id_types.append(pfo_type)
                            mom = pfo.getMomentum()
                            pt = math.sqrt(mom[0]**2 + mom[1]**2)
                            energy = pfo.getEnergy()
                            mis_id_pts.append(pt)
                            mis_id_energies.append(energy)
                            if pfo_type == 11:
                                n_matched_electrons += 1


        if (decay_mode == 4):
            n_electrons.append(n_matched_electrons)
                        
        n_matched_pis_all = n_matched_pis_1p0n + n_matched_unique_pis_3p0n + n_matched_duplicate_pis_3p0n
        if (decay_mode == 0 or decay_mode == 4):
            n_pi_reco.append(n_matched_pis_all)

        if (n_matched_pis_1p0n == 1):
            n_1p0n_reco += 1
            hPiPtMatched1P0N.Fill(vis_pt)
            hPiThetaMatched1P0N.Fill(vis_theta)
            hPiPhiMatched1P0N.Fill(vis_phi)

            if (vis_theta > 0.70 and vis_theta < 2.45):
                hPiPtMatched1P0N_Barrel.Fill(vis_pt)
                hPiPhiMatched1P0N_Barrel.Fill(vis_phi)

            if (vis_theta > 1 and vis_theta < 2):
                hPiPtMatched1P0N_CentBarrel.Fill(vis_pt)
                hPiPhiMatched1P0N_CentBarrel.Fill(vis_phi)

            elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                hPiPtMatched1P0N_Transition.Fill(vis_pt)
                hPiPhiMatched1P0N_Transition.Fill(vis_phi)

            elif (vis_theta < 0.577 or vis_theta > 2.56):
                hPiPtMatched1P0N_Endcap.Fill(vis_pt)
                hPiPhiMatched1P0N_Endcap.Fill(vis_phi)

        elif (n_matched_unique_pis_3p0n == 3 and n_matched_duplicate_pis_3p0n == 0):
            n_3p0n_reco += 1
            hPiPtMatched3P0N.Fill(vis_pt)
            hPiThetaMatched3P0N.Fill(vis_theta)
            hPiPhiMatched3P0N.Fill(vis_phi)

            if (vis_theta > 0.70 and vis_theta < 2.45):
                hPiPtMatched3P0N_Barrel.Fill(vis_pt)
                hPiPhiMatched3P0N_Barrel.Fill(vis_phi)

            if (vis_theta > 1 and vis_theta < 2):
                hPiPtMatched3P0N_CentBarrel.Fill(vis_pt)
                hPiPhiMatched3P0N_CentBarrel.Fill(vis_phi)

            elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                hPiPtMatched3P0N_Transition.Fill(vis_pt)
                hPiPhiMatched3P0N_Transition.Fill(vis_phi)

            elif (vis_theta < 0.577 or vis_theta > 2.56):
                hPiPtMatched3P0N_Endcap.Fill(vis_pt)
                hPiPhiMatched3P0N_Endcap.Fill(vis_phi)
                # if ((65 <= vis_pt and vis_pt <= 128) or (193 <= vis_pt and vis_pt <= 224)):
                    # print(f'\nVisible Theta: {vis_theta}, Visible Pt: {vis_pt}')
                    # for unique_pt in unique_pts:
                        # print(f'Pion Reco Pt: {unique_pt}')

        elif ((n_matched_unique_pis_3p0n != 3 or (n_matched_unique_pis_3p0n == 3 and n_matched_duplicate_pis_3p0n != 0)) and decay_mode == 4):
            if (vis_theta > 1 and vis_theta < 2):
                if (0 <= vis_pt and vis_pt <= 32):
                    true_region = getDetectorRegion(vis_theta)
                    print(f'\nVisible Region: {true_region}, Visible Pt: {vis_pt}')
                    print(f'Event number: {event.getEventNumber()}')
                    print(f'Number of unique reco pis: {n_matched_unique_pis_3p0n}, Number of duplicate reco pis: {n_matched_duplicate_pis_3p0n}')
                    for i in range(3):
                        if true_pis[i] in matched_true_3p0n_pis:
                            reco = 'Yes'
                        else:
                            reco='No'
                        print(f'\nTrue Pion Pt: {true_pi_pts[i]:.2f}, True Pion Energy: {true_pi_energies[i]:.2f}, Reconstructed: {reco}')
                        if (reco == 'No'):
                            daughters = true_pis[i].getDaughters()
                            for daughter in daughters:
                                pdg = abs(daughter.getPDG())
                                notendpoint = daughter.vertexIsNotEndpointOfParent()
                                if (not notendpoint):
                                    pos = daughter.getVertex()
                                    x = pos[0]
                                    y = pos[1]
                                    z = pos[2]
                                    r = math.sqrt(x**2 + y**2 + z**2)
                                    print(f'True Pion Endpoint-Daughter PDG: {pdg}, True Pion Endpoint-Daughter Vertex Radius: {r:.2f}')
                    for k in range(len(unique_pts)):
                        print(f'Pion Reco Pt: {unique_pts[k]}, Pion Reco Energy: {unique_energies[k]}')
                    for l in range(len(duplicate_pts)):
                        print(f'Duplicate Pion Reco Pt: {duplicate_pts[l]}, Duplicate Pion Reco Energy: {duplicate_energies[l]}')
                    for j in range(len(mis_id_types)):
                        print(f'MisID Type: {mis_id_types[j]}, MisID Pt: {mis_id_pts[j]}, MisID Energy: {mis_id_energies[j]}')
            '''
            elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                if (0 <= vis_pt and vis_pt <= 64):
                    true_region = getDetectorRegion(vis_theta)
                    print(f'\nVisible Region: {true_region}, Visible Pt: {vis_pt}')
                    print(f'Number of unique reco pis: {n_matched_unique_pis_3p0n}, Number of duplicate reco pis: {n_matched_duplicate_pis_3p0n}')                    
                    for i in range(3):
                        if true_pis[i] in matched_true_3p0n_pis:
                            reco = 'Yes'
                        else:
                            reco='No'
                        print(f'Pion True Pt: {true_pi_pts[i]}, Pion True Theta: {true_pi_thetas[i]}, Pion True Phi: {true_pi_phis[i]}, Pion True Charge: {true_pi_charges[i]}, Reconstructed: {reco}')
            '''
    reader.close()

print(f'1p0n eff: {n_1p0n_reco/n_1p0n_true}')
print(f'3p0n eff: {n_3p0n_reco/n_3p0n_true}')

# Create confusion matrices
plt.figure()

plt.clf()
disp = ConfusionMatrixDisplay.from_predictions(n_pi_reco, n_pi_true, normalize='pred', include_values=True, cmap='Blues', colorbar=True)
disp.ax_.set_ylabel(r'Number of Reconstructed $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
disp.ax_.set_xlabel(r'Number of True $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
disp.ax_.invert_yaxis()
plt.tight_layout()
plt.savefig('pi_confusion.png')


plt.clf()
pt_bins = np.linspace(0, 320, 11)
bin_indices = np.digitize(tau_pts, pt_bins)
total_electrons_per_bin = []
bin_centers = []
for i in range(1, len(pt_bins)):
    electrons_in_bin = [n_electrons[j] for j in range(len(tau_pts)) if bin_indices[j] == i]
    total = sum(electrons_in_bin)
    total_electrons_per_bin.append(total)
    bin_centers.append(0.5*(pt_bins[i-1]+pt_bins[i]))
plt.bar(bin_centers, total_electrons_per_bin, width=np.diff(pt_bins), align='center', color='m')
plt.xlabel('Transverse Momentum [GeV/c]')
plt.ylabel('Number of Electrons')
plt.tight_layout()
plt.savefig('n_electrons_vs_pt')

# Create 1P0N Pion Eff Hists
hPiPtEff1P0N = hPiPtMatched1P0N.Clone('1p0n_pi_pt_eff')
hPiPtEff1P0N.Divide(hPiPtEff1P0N, hPiPtTrue1P0N, 1, 1, 'B')
hPiPtEff1P0N.SetLineColor(6)
hPiPtEff1P0N.SetLineWidth(2)
hPiPtEff1P0N.SetTitle('1P0N Pion Reconstruction Efficiency vs Pt')
hPiPtEff1P0N.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff1P0N.GetYaxis().SetTitle('#epsilon')
hPiPtEff1P0N.SetStats(0)
hists.append(hPiPtEff1P0N)

hPiThetaEff1P0N = hPiThetaMatched1P0N.Clone('1p0n_pi_theta_eff')
hPiThetaEff1P0N.Divide(hPiThetaEff1P0N, hPiThetaTrue1P0N, 1, 1, 'B')
hPiThetaEff1P0N.SetLineColor(7)
hPiThetaEff1P0N.SetLineWidth(2)
hPiThetaEff1P0N.SetTitle('1P0N Pion Reconstruction Efficiency vs Theta')
hPiThetaEff1P0N.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hPiThetaEff1P0N.GetYaxis().SetTitle('#epsilon')
hPiThetaEff1P0N.SetStats(0)
hists.append(hPiThetaEff1P0N)

hPiPhiEff1P0N = hPiPhiMatched1P0N.Clone('1p0n_pi_phi_eff')
hPiPhiEff1P0N.Divide(hPiPhiEff1P0N, hPiPhiTrue1P0N, 1, 1, 'B')
hPiPhiEff1P0N.SetLineColor(418)
hPiPhiEff1P0N.SetLineWidth(2)
hPiPhiEff1P0N.SetTitle('1P0N Pion Reconstruction Efficiency vs Phi')
hPiPhiEff1P0N.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff1P0N.GetYaxis().SetTitle('#epsilon')
hPiPhiEff1P0N.SetStats(0)
hists.append(hPiPhiEff1P0N)

hPiPtEff1P0N_Barrel = hPiPtMatched1P0N_Barrel.Clone('1p0n_pi_pt_eff_barrel')
hPiPtEff1P0N_Barrel.Divide(hPiPtEff1P0N_Barrel, hPiPtTrue1P0N_Barrel, 1, 1, 'B')
hPiPtEff1P0N_Barrel.SetLineColor(6)
hPiPtEff1P0N_Barrel.SetLineWidth(2)
hPiPtEff1P0N_Barrel.SetTitle('1P0N Pion Reconstruction Efficiency vs Pt (Barrel)')
hPiPtEff1P0N_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff1P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPiPtEff1P0N_Barrel.SetStats(0)
hists.append(hPiPtEff1P0N_Barrel)

hPiPhiEff1P0N_Barrel = hPiPhiMatched1P0N_Barrel.Clone('1p0n_pi_phi_eff_barrel')
hPiPhiEff1P0N_Barrel.Divide(hPiPhiEff1P0N_Barrel, hPiPhiTrue1P0N_Barrel, 1, 1, 'B')
hPiPhiEff1P0N_Barrel.SetLineColor(418)
hPiPhiEff1P0N_Barrel.SetLineWidth(2)
hPiPhiEff1P0N_Barrel.SetTitle('1P0N Pion Reconstruction Efficiency vs Phi (Barrel)')
hPiPhiEff1P0N_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff1P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPiPhiEff1P0N_Barrel.SetStats(0)
hists.append(hPiPhiEff1P0N_Barrel)

hPiPtEff1P0N_CentBarrel = hPiPtMatched1P0N_CentBarrel.Clone('1p0n_pi_pt_eff_centbarrel')
hPiPtEff1P0N_CentBarrel.Divide(hPiPtEff1P0N_CentBarrel, hPiPtTrue1P0N_CentBarrel, 1, 1, 'B')
hPiPtEff1P0N_CentBarrel.SetLineColor(6)
hPiPtEff1P0N_CentBarrel.SetLineWidth(2)
hPiPtEff1P0N_CentBarrel.SetTitle('1P0N Pion Reconstruction Efficiency vs Pt (Central Barrel)')
hPiPtEff1P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff1P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPiPtEff1P0N_CentBarrel.SetStats(0)
hists.append(hPiPtEff1P0N_CentBarrel)

hPiPhiEff1P0N_CentBarrel = hPiPhiMatched1P0N_CentBarrel.Clone('1p0n_pi_phi_eff_centbarrel')
hPiPhiEff1P0N_CentBarrel.Divide(hPiPhiEff1P0N_CentBarrel, hPiPhiTrue1P0N_CentBarrel, 1, 1, 'B')
hPiPhiEff1P0N_CentBarrel.SetLineColor(418)
hPiPhiEff1P0N_CentBarrel.SetLineWidth(2)
hPiPhiEff1P0N_CentBarrel.SetTitle('1P0N Pion Reconstruction Efficiency vs Phi (Central Barrel)')
hPiPhiEff1P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff1P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPiPhiEff1P0N_CentBarrel.SetStats(0)
hists.append(hPiPhiEff1P0N_CentBarrel)

hPiPtEff1P0N_Endcap = hPiPtMatched1P0N_Endcap.Clone('1p0n_pi_pt_eff_endcap')
hPiPtEff1P0N_Endcap.Divide(hPiPtEff1P0N_Endcap, hPiPtTrue1P0N_Endcap, 1, 1, 'B')
hPiPtEff1P0N_Endcap.SetLineColor(6)
hPiPtEff1P0N_Endcap.SetLineWidth(2)
hPiPtEff1P0N_Endcap.SetTitle('1P0N Pion Reconstruction Efficiency vs Pt (Endcap)')
hPiPtEff1P0N_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff1P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPiPtEff1P0N_Endcap.SetStats(0)
hists.append(hPiPtEff1P0N_Endcap)

hPiPhiEff1P0N_Endcap = hPiPhiMatched1P0N_Endcap.Clone('1p0n_pi_phi_eff_endcap')
hPiPhiEff1P0N_Endcap.Divide(hPiPhiEff1P0N_Endcap, hPiPhiTrue1P0N_Endcap, 1, 1, 'B')
hPiPhiEff1P0N_Endcap.SetLineColor(418)
hPiPhiEff1P0N_Endcap.SetLineWidth(2)
hPiPhiEff1P0N_Endcap.SetTitle('1P0N Pion Reconstruction Efficiency vs Phi (Endcap)')
hPiPhiEff1P0N_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff1P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPiPhiEff1P0N_Endcap.SetStats(0)
hists.append(hPiPhiEff1P0N_Endcap)

hPiPtEff1P0N_Transition = hPiPtMatched1P0N_Transition.Clone('1p0n_pi_pt_eff_transition')
hPiPtEff1P0N_Transition.Divide(hPiPtEff1P0N_Transition, hPiPtTrue1P0N_Transition, 1, 1, 'B')
hPiPtEff1P0N_Transition.SetLineColor(6)
hPiPtEff1P0N_Transition.SetLineWidth(2)
hPiPtEff1P0N_Transition.SetTitle('1P0N Pion Reconstruction Efficiency vs Pt (Transition)')
hPiPtEff1P0N_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff1P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPiPtEff1P0N_Transition.SetStats(0)
hists.append(hPiPtEff1P0N_Transition)

hPiPhiEff1P0N_Transition = hPiPhiMatched1P0N_Transition.Clone('1p0n_pi_phi_eff_transition')
hPiPhiEff1P0N_Transition.Divide(hPiPhiEff1P0N_Transition, hPiPhiTrue1P0N_Transition, 1, 1, 'B')
hPiPhiEff1P0N_Transition.SetLineColor(418)
hPiPhiEff1P0N_Transition.SetLineWidth(2)
hPiPhiEff1P0N_Transition.SetTitle('1P0N Pion Reconstruction Efficiency vs Phi (Transition)')
hPiPhiEff1P0N_Transition.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff1P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPiPhiEff1P0N_Transition.SetStats(0)
hists.append(hPiPhiEff1P0N_Transition)

# Create 3P0N Pion Eff Hists
hPiPtEff3P0N = hPiPtMatched3P0N.Clone('3p0n_pi_pt_eff')
hPiPtEff3P0N.Divide(hPiPtEff3P0N, hPiPtTrue3P0N, 1, 1, 'B')
hPiPtEff3P0N.SetLineColor(6)
hPiPtEff3P0N.SetLineWidth(2)
hPiPtEff3P0N.SetTitle('3P0N Pion Reconstruction Efficiency vs Pt')
hPiPtEff3P0N.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff3P0N.GetYaxis().SetTitle('#epsilon')
hPiPtEff3P0N.SetStats(0)
hists.append(hPiPtEff3P0N)

hPiThetaEff3P0N = hPiThetaMatched3P0N.Clone('3p0n_pi_theta_eff')
hPiThetaEff3P0N.Divide(hPiThetaEff3P0N, hPiThetaTrue3P0N, 1, 1, 'B')
hPiThetaEff3P0N.SetLineColor(7)
hPiThetaEff3P0N.SetLineWidth(2)
hPiThetaEff3P0N.SetTitle('3P0N Pion Reconstruction Efficiency vs Theta')
hPiThetaEff3P0N.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hPiThetaEff3P0N.GetYaxis().SetTitle('#epsilon')
hPiThetaEff3P0N.SetStats(0)
hists.append(hPiThetaEff3P0N)

hPiPhiEff3P0N = hPiPhiMatched3P0N.Clone('3p0n_pi_phi_eff')
hPiPhiEff3P0N.Divide(hPiPhiEff3P0N, hPiPhiTrue3P0N, 1, 1, 'B')
hPiPhiEff3P0N.SetLineColor(438)
hPiPhiEff3P0N.SetLineWidth(2)
hPiPhiEff3P0N.SetTitle('3P0N Pion Reconstruction Efficiency vs Phi')
hPiPhiEff3P0N.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff3P0N.GetYaxis().SetTitle('#epsilon')
hPiPhiEff3P0N.SetStats(0)
hists.append(hPiPhiEff3P0N)

hPiPtEff3P0N_Barrel = hPiPtMatched3P0N_Barrel.Clone('3p0n_pi_pt_eff_barrel')
hPiPtEff3P0N_Barrel.Divide(hPiPtEff3P0N_Barrel, hPiPtTrue3P0N_Barrel, 1, 1, 'B')
hPiPtEff3P0N_Barrel.SetLineColor(6)
hPiPtEff3P0N_Barrel.SetLineWidth(2)
hPiPtEff3P0N_Barrel.SetTitle('3P0N Pion Reconstruction Efficiency vs Pt (Barrel)')
hPiPtEff3P0N_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff3P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPiPtEff3P0N_Barrel.SetStats(0)
hists.append(hPiPtEff3P0N_Barrel)

hPiPhiEff3P0N_Barrel = hPiPhiMatched3P0N_Barrel.Clone('3p0n_pi_phi_eff_barrel')
hPiPhiEff3P0N_Barrel.Divide(hPiPhiEff3P0N_Barrel, hPiPhiTrue3P0N_Barrel, 1, 1, 'B')
hPiPhiEff3P0N_Barrel.SetLineColor(418)
hPiPhiEff3P0N_Barrel.SetLineWidth(2)
hPiPhiEff3P0N_Barrel.SetTitle('3P0N Pion Reconstruction Efficiency vs Phi (Barrel)')
hPiPhiEff3P0N_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff3P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPiPhiEff3P0N_Barrel.SetStats(0)
hists.append(hPiPhiEff3P0N_Barrel)

hPiPtEff3P0N_CentBarrel = hPiPtMatched3P0N_CentBarrel.Clone('3p0n_pi_pt_eff_centbarrel')
hPiPtEff3P0N_CentBarrel.Divide(hPiPtEff3P0N_CentBarrel, hPiPtTrue3P0N_CentBarrel, 1, 1, 'B')
hPiPtEff3P0N_CentBarrel.SetLineColor(6)
hPiPtEff3P0N_CentBarrel.SetLineWidth(2)
hPiPtEff3P0N_CentBarrel.SetTitle('3P0N Pion Reconstruction Efficiency vs Pt (Central Barrel)')
hPiPtEff3P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff3P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPiPtEff3P0N_CentBarrel.SetStats(0)
hists.append(hPiPtEff3P0N_CentBarrel)

hPiPhiEff3P0N_CentBarrel = hPiPhiMatched3P0N_CentBarrel.Clone('3p0n_pi_phi_eff_centbarrel')
hPiPhiEff3P0N_CentBarrel.Divide(hPiPhiEff3P0N_CentBarrel, hPiPhiTrue3P0N_CentBarrel, 1, 1, 'B')
hPiPhiEff3P0N_CentBarrel.SetLineColor(418)
hPiPhiEff3P0N_CentBarrel.SetLineWidth(2)
hPiPhiEff3P0N_CentBarrel.SetTitle('3P0N Pion Reconstruction Efficiency vs Phi (Central Barrel)')
hPiPhiEff3P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff3P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPiPhiEff3P0N_CentBarrel.SetStats(0)
hists.append(hPiPhiEff3P0N_CentBarrel)

hPiPtEff3P0N_Endcap = hPiPtMatched3P0N_Endcap.Clone('3p0n_pi_pt_eff_endcap')
hPiPtEff3P0N_Endcap.Divide(hPiPtEff3P0N_Endcap, hPiPtTrue3P0N_Endcap, 1, 1, 'B')
hPiPtEff3P0N_Endcap.SetLineColor(6)
hPiPtEff3P0N_Endcap.SetLineWidth(2)
hPiPtEff3P0N_Endcap.SetTitle('3P0N Pion Reconstruction Efficiency vs Pt (Endcap)')
hPiPtEff3P0N_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff3P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPiPtEff3P0N_Endcap.SetStats(0)
hists.append(hPiPtEff3P0N_Endcap)

hPiPhiEff3P0N_Endcap = hPiPhiMatched3P0N_Endcap.Clone('3p0n_pi_phi_eff_endcap')
hPiPhiEff3P0N_Endcap.Divide(hPiPhiEff3P0N_Endcap, hPiPhiTrue3P0N_Endcap, 1, 1, 'B')
hPiPhiEff3P0N_Endcap.SetLineColor(418)
hPiPhiEff3P0N_Endcap.SetLineWidth(2)
hPiPhiEff3P0N_Endcap.SetTitle('3P0N Pion Reconstruction Efficiency vs Phi (Endcap)')
hPiPhiEff3P0N_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff3P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPiPhiEff3P0N_Endcap.SetStats(0)
hists.append(hPiPhiEff3P0N_Endcap)

hPiPtEff3P0N_Transition = hPiPtMatched3P0N_Transition.Clone('3p0n_pi_pt_eff_transition')
hPiPtEff3P0N_Transition.Divide(hPiPtEff3P0N_Transition, hPiPtTrue3P0N_Transition, 1, 1, 'B')
hPiPtEff3P0N_Transition.SetLineColor(6)
hPiPtEff3P0N_Transition.SetLineWidth(2)
hPiPtEff3P0N_Transition.SetTitle('3P0N Pion Reconstruction Efficiency vs Pt (Transition)')
hPiPtEff3P0N_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPiPtEff3P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPiPtEff3P0N_Transition.SetStats(0)
hists.append(hPiPtEff3P0N_Transition)

hPiPhiEff3P0N_Transition = hPiPhiMatched3P0N_Transition.Clone('3p0n_pi_phi_eff_transition')
hPiPhiEff3P0N_Transition.Divide(hPiPhiEff3P0N_Transition, hPiPhiTrue3P0N_Transition, 1, 1, 'B')
hPiPhiEff3P0N_Transition.SetLineColor(418)
hPiPhiEff3P0N_Transition.SetLineWidth(2)
hPiPhiEff3P0N_Transition.SetTitle('3P0N Pion Reconstruction Efficiency vs Phi (Transition)')
hPiPhiEff3P0N_Transition.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPiPhiEff3P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPiPhiEff3P0N_Transition.SetStats(0)
hists.append(hPiPhiEff3P0N_Transition)

# Write to output file
output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
    hist.Write()
output_file.Close()

# Draw hists and save as PNG
for hist in hists:
    filename = hist.GetName() + '.png'
    canvas = TCanvas()
    hist.Draw()
    canvas.SaveAs(filename)
