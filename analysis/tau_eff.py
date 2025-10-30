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

# Get tau decay mode
def getDecayMode(mcTau):

    '''
    Decay modes for tau-:
    0: -211, 16, (pi-, nu-tau) (10.82%)
    1: -211, 16, 111 (pi-, pi0, nu-tau) (25.49%)
    2: -211, 16, 111, 111 (pi-, pi0x2, nu-tau) (9.26%)
    3: -211, 16, 111, 111, 111 (pi-, pi0x3, nu-tau) (1.04%)
    4: -211, -211, 16, 211 (3-prong plus nu-tau) (8.99%)
    5: -211, -211, 16, 111, 211 (3-prong plus pi0 and nu-tau) (2.74%)
    6: -12, 11, 16 (nu-e-bar, e-, nu-tau) (17.82%)
    7: -14, 13, 16 (nu-mu-bar, mu-, nu-tau) (17.39%) 
    '''
    # Initialize list to store daughter pdgs
    daughter_pdgs = []

    n_daughters = len(mcTau.getDaughters())
    
    # Loop over daughters and store pdgs
    for daughter in mcTau.getDaughters():
        daughter_pdgs.append(daughter.getPDG())

    if (mcTau.getPDG() == 15):

        # Sort daughter pdgs from lowest to highest
        daughter_pdgs.sort()

        # Hadronic decays
        if (daughter_pdgs[0] == -211):
            if (daughter_pdgs[1] == -211):
                # 3-prong
                if (daughter_pdgs[3] == 211):
                    # 3-prong no neutrals
                    return 4
                elif (daughter_pdgs[3] == 111):
                    # 3-prong with neutral
                    return 5
            elif (daughter_pdgs[1] == 16):
                # 1-prong
                if (n_daughters == 2):
                    # No neutrals
                    return 0
                elif (daughter_pdgs[2] == 111):
                    if (n_daughters == 3):
                        # 1 neutral
                        return 1
                    elif (n_daughters == 4):
                        # 2 neutrals
                        return 2
                    elif (n_daughters == 5):
                        # 3 neutrals
                        return 3
        # Leptonic decays
        elif (daughter_pdgs[0] == -12):
            # Electron
            return 6
        else:
            # Muon
            return 7

# Get visible 4-momentum
def getVisibleP(mc_tau):

    E_vis = 0
    px_vis = 0
    py_vis = 0
    pz_vis = 0
    
    # Loop through mc daughters
    mc_daughters = mc_tau.getDaughters()
    for mc_daughter in mc_daughters:

        # Ignore neutrinos
        pdg = abs(mc_daughter.getPDG())
        if (pdg == 12 or pdg == 14 or pdg == 16):
            continue

        # Sum visible observables
        E_vis += mc_daughter.getEnergy()
        px_vis += mc_daughter.getMomentum()[0]
        py_vis += mc_daughter.getMomentum()[1]
        pz_vis += mc_daughter.getMomentum()[2]

    # 4-momentum vector
    vis_properties = [E_vis, px_vis, py_vis, pz_vis]
    return vis_properties

# Get transverse momentum
def getPt(vis_props):
    return math.sqrt(vis_props[1]**2 + vis_props[2]**2)

# Get 3-momentum magnitude
def getP(vis_props):
    return math.sqrt(vis_props[1]**2 + vis_props[2]**2 + vis_props[3]**2)

# Get polar angle
def getTheta(vis_props):
    p = getP(vis_props)
    pz = vis_props[3]
    return math.acos(pz/p)

# Get azimuthal angle
def getPhi(vis_props):
    px = vis_props[1]
    pt = getPt(vis_props)
    return math.acos(px/pt)

# Check if particle is a direct product of tau decay
def isTauDaughter(mc_pi):

    # Loop through parents
    parents = mc_pi.getParents()
    for parent in parents:
        # Check if any parent is a tau
        if (abs(parent.getPDG()) == 15):
            return True
    return False

# Get number of matched reconstructed charged pions in a reconstructed tau
def getNRecoChargedPis(reco_tau, mc_particles):
    matched_mc_pis = []
    
    # Loop through reco tau pfos
    pfos = reco_tau.getParticles()
    for pfo in pfos:

        # Tag charged pions
        if (abs(pfo.getType()) == 211):

            # Get observables
            pfo_mom = pfo.getMomentum()
            pfo_px = pfo_mom[0]
            pfo_py = pfo_mom[1]
            pfo_pz = pfo_mom[2]
            pfo_pt = math.sqrt(pfo_px**2 + pfo_py**2)
            pfo_theta = math.atan2(pfo_pt, pfo_pz)
            # pfo_eta = -math.log(math.tan(theta/2))
            pfo_phi = math.atan2(pfo_py, pfo_px)
            pfo_energy = pfo.getEnergy()

            # Loop over mc particles
            deltaR_mc_pis = []
            for mc_particle in mc_particles:

                # Tag charged pions from tau decay
                if (abs(mc_particle.getPDG()) == 211 and isTauDaughter(mc_particle)):
                            
                    # Get observables
                    mc_mom = mc_particle.getMomentum()
                    mc_px = mc_mom[0]
                    mc_py = mc_mom[1]
                    mc_pz = mc_mom[2]
                    mc_pt = math.sqrt(mc_px**2 + mc_py**2)
                    mc_theta = math.atan2(mc_pt, mc_pz)
                    # mc_eta = -math.log(math.tan(theta/2))
                    mc_phi = math.atan2(mc_py, mc_px)
                    
                    # Compute deltaR
                    deltaTheta = pfo_theta - mc_theta
                    deltaPhi = pfo_phi - mc_phi
                    while (deltaPhi > math.pi):
                        deltaPhi -= 2*math.pi
                    while (deltaPhi < -math.pi):
                        deltaPhi += 2*math.pi
                    deltaR = math.sqrt(deltaTheta**2 + deltaPhi**2)
                    
                    # Save mc particle if deltaR is less than 0.1 rad
                    if (deltaR < 0.1):
                        deltaR_mc_pis.append(mc_particle)
                        
            # Loop over mc particles with deltaR less than 0.1 rad
            min_energy_dif = 1e6
            matched_mc_pi = None
            for deltaR_mc_pi in deltaR_mc_pis:
                
                # Compute energy difference between mc particle and pfo
                mc_energy = deltaR_mc_pi.getEnergy()
                energy_dif = abs(pfo_energy - mc_energy)

                # MC particle matched to pfo if it has smallest energy difference
                if (energy_dif < min_energy_dif):
                    min_energy_dif = energy_dif
                    matched_mc_pi = deltaR_mc_pi

            if (matched_mc_pi is not None):
                matched_mc_pis.append(matched_mc_pi)

    return len(matched_mc_pis)

# Get number of reconstructed neutral pions in a reconstructed tau
def getNRecoNeutralPis(reco_tau):

    # Count number of neutral pions by counting number of photons
    
    n_reco_pi0s = 0

    # Loop through reco tau pfos
    pfos = reco_tau.getParticles()
    for pfo in pfos:

        # Tag photons
        pfo_type = abs(pfo.getType())
        if (pfo_type == 22):
            n_reco_pi0s += 1

    return n_reco_pi0s

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='tau_eff.root')

# Transverse momentum lower bound
parser.add_argument('--pTMin', type=float, default=0)

# Transverse momentum upper bound
parser.add_argument('--pTMax', type=float, default=250)

# Polar angle lower bound
parser.add_argument('--thetaMin', type=float, default=0)

# Polar angle upper bound
parser.add_argument('--thetaMax', type=float, default=math.pi)

# Azimuthal angle lower bound
parser.add_argument('--phiMin', type=float, default=0)

# Azimuthal angle upper bound
parser.add_argument('--phiMax', type=float, default=math.pi)

args = parser.parse_args()

# Initialize histograms
hists = []

# 1P0N Matched
hMatchedTauPt1P0N = TH1F('matched_1p0n_tau_pt', 'Matched 1P0N Tau Pt', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P0N)

hMatchedTauTheta1P0N = TH1F('matched_1p0n_tau_theta', 'Matched 1P0N Tau Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hMatchedTauTheta1P0N)

hMatchedTauPhi1P0N = TH1F('matched_1p0n_tau_phi', 'Matched 1P0N Tau Phi', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P0N)

hMatchedTauPt1P0N_Barrel = TH1F('matched_1p0n_tau_pt_barrel', 'Matched 1P0N Tau Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P0N_Barrel)

hMatchedTauPhi1P0N_Barrel = TH1F('matched_1p0n_tau_phi_barrel', 'Matched 1P0N Tau Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P0N_Barrel)

hMatchedTauPt1P0N_CentBarrel = TH1F('matched_1p0n_tau_pt_centbarrel', 'Matched 1P0N Tau Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P0N_CentBarrel)

hMatchedTauPhi1P0N_CentBarrel = TH1F('matched_1p0n_tau_phi_centbarrel', 'Matched 1P0N Tau Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P0N_CentBarrel)

hMatchedTauPt1P0N_Endcap = TH1F('matched_1p0n_tau_pt_endcap', 'Matched 1P0N Tau Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P0N_Endcap)

hMatchedTauPhi1P0N_Endcap = TH1F('matched_1p0n_tau_phi_endcap', 'Matched 1P0N Tau Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P0N_Endcap)

hMatchedTauPt1P0N_Transition = TH1F('matched_1p0n_tau_pt_transition', 'Matched 1P0N Tau Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P0N_Transition)

hMatchedTauPhi1P0N_Transition = TH1F('matched_1p0n_tau_phi_transition', 'Matched 1P0N Tau Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P0N_Transition)

# 1P1N Matched
hMatchedTauPt1P1N = TH1F('matched_1p1n_tau_pt', 'Matched 1P1N Tau Pt', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P1N)

hMatchedTauTheta1P1N = TH1F('matched_1p1n_tau_theta', 'Matched 1P1N Tau Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hMatchedTauTheta1P1N)

hMatchedTauPhi1P1N = TH1F('matched_1p1n_tau_phi', 'Matched 1P1N Tau Phi', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P1N)

hMatchedTauPt1P1N_Barrel = TH1F('matched_1p1n_tau_pt_barrel', 'Matched 1P1N Tau Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P1N_Barrel)

hMatchedTauPhi1P1N_Barrel = TH1F('matched_1p1n_tau_phi_barrel', 'Matched 1P1N Tau Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P1N_Barrel)

hMatchedTauPt1P1N_CentBarrel = TH1F('matched_1p1n_tau_pt_centbarrel', 'Matched 1P1N Tau Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P1N_CentBarrel)

hMatchedTauPhi1P1N_CentBarrel = TH1F('matched_1p1n_tau_phi_centbarrel', 'Matched 1P1N Tau Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P1N_CentBarrel)

hMatchedTauPt1P1N_Endcap = TH1F('matched_1p1n_tau_pt_endcap', 'Matched 1P1N Tau Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P1N_Endcap)

hMatchedTauPhi1P1N_Endcap = TH1F('matched_1p1n_tau_phi_endcap', 'Matched 1P1N Tau Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P1N_Endcap)

hMatchedTauPt1P1N_Transition = TH1F('matched_1p1n_tau_pt_transition', 'Matched 1P1N Tau Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P1N_Transition)

hMatchedTauPhi1P1N_Transition = TH1F('matched_1p1n_tau_phi_transition', 'Matched 1P1N Tau Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P1N_Transition)

# 1P2N Matched
hMatchedTauPt1P2N = TH1F('matched_1p2n_tau_pt', 'Matched 1P2N Tau Pt', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P2N)

hMatchedTauTheta1P2N = TH1F('matched_1p2n_tau_theta', 'Matched 1P2N Tau Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hMatchedTauTheta1P2N)

hMatchedTauPhi1P2N = TH1F('matched_1p2n_tau_phi', 'Matched 1P2N Tau Phi', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P2N)

hMatchedTauPt1P2N_Barrel = TH1F('matched_1p2n_tau_pt_barrel', 'Matched 1P2N Tau Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P2N_Barrel)

hMatchedTauPhi1P2N_Barrel = TH1F('matched_1p2n_tau_phi_barrel', 'Matched 1P2N Tau Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P2N_Barrel)

hMatchedTauPt1P2N_CentBarrel = TH1F('matched_1p2n_tau_pt_centbarrel', 'Matched 1P2N Tau Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P2N_CentBarrel)

hMatchedTauPhi1P2N_CentBarrel = TH1F('matched_1p2n_tau_phi_centbarrel', 'Matched 1P2N Tau Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P2N_CentBarrel)

hMatchedTauPt1P2N_Endcap = TH1F('matched_1p2n_tau_pt_endcap', 'Matched 1P2N Tau Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P2N_Endcap)

hMatchedTauPhi1P2N_Endcap = TH1F('matched_1p2n_tau_phi_endcap', 'Matched 1P2N Tau Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P2N_Endcap)

hMatchedTauPt1P2N_Transition = TH1F('matched_1p2n_tau_pt_transition', 'Matched 1P2N Tau Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1P2N_Transition)

hMatchedTauPhi1P2N_Transition = TH1F('matched_1p2n_tau_phi_transition', 'Matched 1P2N Tau Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1P2N_Transition)

# 1PXN Matched
hMatchedTauPt1PXN = TH1F('matched_1pXn_tau_pt', 'Matched 1PXN Tau Pt', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1PXN)

hMatchedTauTheta1PXN = TH1F('matched_1pXn_tau_theta', 'Matched 1PXN Tau Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hMatchedTauTheta1PXN)

hMatchedTauPhi1PXN = TH1F('matched_1pXn_tau_phi', 'Matched 1PXN Tau Phi', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1PXN)

hMatchedTauPt1PXN_Barrel = TH1F('matched_1pXn_tau_pt_barrel', 'Matched 1PXN Tau Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1PXN_Barrel)

hMatchedTauPhi1PXN_Barrel = TH1F('matched_1pXn_tau_phi_barrel', 'Matched 1PXN Tau Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1PXN_Barrel)

hMatchedTauPt1PXN_CentBarrel = TH1F('matched_1pXn_tau_pt_centbarrel', 'Matched 1PXN Tau Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1PXN_CentBarrel)

hMatchedTauPhi1PXN_CentBarrel = TH1F('matched_1pXn_tau_phi_centbarrel', 'Matched 1PXN Tau Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1PXN_CentBarrel)

hMatchedTauPt1PXN_Endcap = TH1F('matched_1pXn_tau_pt_endcap', 'Matched 1PXN Tau Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1PXN_Endcap)

hMatchedTauPhi1PXN_Endcap = TH1F('matched_1pXn_tau_phi_endcap', 'Matched 1PXN Tau Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1PXN_Endcap)

hMatchedTauPt1PXN_Transition = TH1F('matched_1pXn_tau_pt_transition', 'Matched 1PXN Tau Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt1PXN_Transition)

hMatchedTauPhi1PXN_Transition = TH1F('matched_1pXn_tau_phi_transition', 'Matched 1PXN Tau Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi1PXN_Transition)

# 3P0N Matched
hMatchedTauPt3P0N = TH1F('matched_3p0n_tau_pt', 'Matched 3P0N Tau Pt', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt3P0N)

hMatchedTauTheta3P0N = TH1F('matched_3p0n_tau_theta', 'Matched 3P0N Tau Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hMatchedTauTheta3P0N)

hMatchedTauPhi3P0N = TH1F('matched_3p0n_tau_phi', 'Matched 3P0N Tau Phi', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi3P0N)

hMatchedTauPt3P0N_Barrel = TH1F('matched_3p0n_tau_pt_barrel', 'Matched 3P0N Tau Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt3P0N_Barrel)

hMatchedTauPhi3P0N_Barrel = TH1F('matched_3p0n_tau_phi_barrel', 'Matched 3P0N Tau Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi3P0N_Barrel)

hMatchedTauPt3P0N_CentBarrel = TH1F('matched_3p0n_tau_pt_centbarrel', 'Matched 3P0N Tau Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt3P0N_CentBarrel)

hMatchedTauPhi3P0N_CentBarrel = TH1F('matched_3p0n_tau_phi_centbarrel', 'Matched 3P0N Tau Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi3P0N_CentBarrel)

hMatchedTauPt3P0N_Endcap = TH1F('matched_3p0n_tau_pt_endcap', 'Matched 3P0N Tau Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt3P0N_Endcap)

hMatchedTauPhi3P0N_Endcap = TH1F('matched_3p0n_tau_phi_endcap', 'Matched 3P0N Tau Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi3P0N_Endcap)

hMatchedTauPt3P0N_Transition = TH1F('matched_3p0n_tau_pt_transition', 'Matched 3P0N Tau Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hMatchedTauPt3P0N_Transition)

hMatchedTauPhi3P0N_Transition = TH1F('matched_3p0n_tau_phi_transition', 'Matched 3P0N Tau Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hMatchedTauPhi3P0N_Transition)

# 1P0N True
hTrueTauVisPt1P0N = TH1F('true_1p0n_tau_vis_pt', 'True 1P0N Tau Visible Pt', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P0N)

hTrueTauVisTheta1P0N = TH1F('true_1p0n_tau_vis_theta', 'True 1P0N Tau Visible Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hTrueTauVisTheta1P0N)

hTrueTauVisPhi1P0N = TH1F('true_1p0n_tau_vis_phi', 'True 1P0N Tau Visible Phi', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P0N)

hTrueTauVisPt1P0N_Barrel = TH1F('true_1p0n_tau_vis_pt_barrel', 'True 1P0N Tau Visible Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P0N_Barrel)

hTrueTauVisPhi1P0N_Barrel = TH1F('true_1p0n_tau_vis_phi_barrel', 'True 1P0N Tau Visible Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P0N_Barrel)

hTrueTauVisPt1P0N_CentBarrel = TH1F('true_1p0n_tau_vis_pt_centbarrel', 'True 1P0N Tau Visible Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P0N_CentBarrel)

hTrueTauVisPhi1P0N_CentBarrel = TH1F('true_1p0n_tau_vis_phi_centbarrel', 'True 1P0N Tau Visible Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P0N_CentBarrel)

hTrueTauVisPt1P0N_Endcap = TH1F('true_1p0n_tau_vis_pt_endcap', 'True 1P0N Tau Visible Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P0N_Endcap)

hTrueTauVisPhi1P0N_Endcap = TH1F('true_1p0n_tau_vis_phi_endcap', 'True 1P0N Tau Visible Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P0N_Endcap)

hTrueTauVisPt1P0N_Transition = TH1F('true_1p0n_tau_vis_pt_transition', 'True 1P0N Tau Visible Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P0N_Transition)

hTrueTauVisPhi1P0N_Transition = TH1F('true_1p0n_tau_vis_phi_transition', 'True 1P0N Tau Visible Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P0N_Transition)

# 1P1N True
hTrueTauVisPt1P1N = TH1F('true_1p1n_tau_vis_pt', 'True 1P1N Tau Visible Pt', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P1N)

hTrueTauVisTheta1P1N = TH1F('true_1p1n_tau_vis_theta', 'True 1P1N Tau Visible Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hTrueTauVisTheta1P1N)

hTrueTauVisPhi1P1N = TH1F('true_1p1n_tau_vis_phi', 'True 1P1N Tau Visible Phi', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P1N)

hTrueTauVisPt1P1N_Barrel = TH1F('true_1p1n_tau_vis_pt_barrel', 'True 1P1N Tau Visible Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P1N_Barrel)

hTrueTauVisPhi1P1N_Barrel = TH1F('true_1p1n_tau_vis_phi_barrel', 'True 1P1N Tau Visible Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P1N_Barrel)

hTrueTauVisPt1P1N_CentBarrel = TH1F('true_1p1n_tau_vis_pt_centbarrel', 'True 1P1N Tau Visible Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P1N_CentBarrel)

hTrueTauVisPhi1P1N_CentBarrel = TH1F('true_1p1n_tau_vis_phi_centbarrel', 'True 1P1N Tau Visible Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P1N_CentBarrel)

hTrueTauVisPt1P1N_Endcap = TH1F('true_1p1n_tau_vis_pt_endcap', 'True 1P1N Tau Visible Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P1N_Endcap)

hTrueTauVisPhi1P1N_Endcap = TH1F('true_1p1n_tau_vis_phi_endcap', 'True 1P1N Tau Visible Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P1N_Endcap)

hTrueTauVisPt1P1N_Transition = TH1F('true_1p1n_tau_vis_pt_transition', 'True 1P1N Tau Visible Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P1N_Transition)

hTrueTauVisPhi1P1N_Transition = TH1F('true_1p1n_tau_vis_phi_transition', 'True 1P1N Tau Visible Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P1N_Transition)

# 1P2N True
hTrueTauVisPt1P2N = TH1F('true_1p2n_tau_vis_pt', 'True 1P2N Tau Visible Pt', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P2N)

hTrueTauVisTheta1P2N = TH1F('true_1p2n_tau_vis_theta', 'True 1P2N Tau Visible Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hTrueTauVisTheta1P2N)

hTrueTauVisPhi1P2N = TH1F('true_1p2n_tau_vis_phi', 'True 1P2N Tau Visible Phi', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P2N)

hTrueTauVisPt1P2N_Barrel = TH1F('true_1p2n_tau_vis_pt_barrel', 'True 1P2N Tau Visible Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P2N_Barrel)

hTrueTauVisPhi1P2N_Barrel = TH1F('true_1p2n_tau_vis_phi_barrel', 'True 1P2N Tau Visible Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P2N_Barrel)

hTrueTauVisPt1P2N_CentBarrel = TH1F('true_1p2n_tau_vis_pt_centbarrel', 'True 1P2N Tau Visible Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P2N_CentBarrel)

hTrueTauVisPhi1P2N_CentBarrel = TH1F('true_1p2n_tau_vis_phi_centbarrel', 'True 1P2N Tau Visible Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P2N_CentBarrel)

hTrueTauVisPt1P2N_Endcap = TH1F('true_1p2n_tau_vis_pt_endcap', 'True 1P2N Tau Visible Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P2N_Endcap)

hTrueTauVisPhi1P2N_Endcap = TH1F('true_1p2n_tau_vis_phi_endcap', 'True 1P2N Tau Visible Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P2N_Endcap)

hTrueTauVisPt1P2N_Transition = TH1F('true_1p2n_tau_vis_pt_transition', 'True 1P2N Tau Visible Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1P2N_Transition)

hTrueTauVisPhi1P2N_Transition = TH1F('true_1p2n_tau_vis_phi_transition', 'True 1P2N Tau Visible Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1P2N_Transition)

# 1PXN True
hTrueTauVisPt1PXN = TH1F('true_1pXn_tau_vis_pt', 'True 1PXN Tau Visible Pt', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1PXN)

hTrueTauVisTheta1PXN = TH1F('true_1pXn_tau_vis_theta', 'True 1PXN Tau Visible Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hTrueTauVisTheta1PXN)

hTrueTauVisPhi1PXN = TH1F('true_1pXn_tau_vis_phi', 'True 1PXN Tau Visible Phi', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1PXN)

hTrueTauVisPt1PXN_Barrel = TH1F('true_1pXn_tau_vis_pt_barrel', 'True 1PXN Tau Visible Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1PXN_Barrel)

hTrueTauVisPhi1PXN_Barrel = TH1F('true_1pXn_tau_vis_phi_barrel', 'True 1PXN Tau Visible Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1PXN_Barrel)

hTrueTauVisPt1PXN_CentBarrel = TH1F('true_1pXn_tau_vis_pt_centbarrel', 'True 1PXN Tau Visible Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1PXN_CentBarrel)

hTrueTauVisPhi1PXN_CentBarrel = TH1F('true_1pXn_tau_vis_phi_centbarrel', 'True 1PXN Tau Visible Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1PXN_CentBarrel)

hTrueTauVisPt1PXN_Endcap = TH1F('true_1pXn_tau_vis_pt_endcap', 'True 1PXN Tau Visible Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1PXN_Endcap)

hTrueTauVisPhi1PXN_Endcap = TH1F('true_1pXn_tau_vis_phi_endcap', 'True 1PXN Tau Visible Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1PXN_Endcap)

hTrueTauVisPt1PXN_Transition = TH1F('true_1pXn_tau_vis_pt_transition', 'True 1PXN Tau Visible Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt1PXN_Transition)

hTrueTauVisPhi1PXN_Transition = TH1F('true_1pXn_tau_vis_phi_transition', 'True 1PXN Tau Visible Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi1PXN_Transition)

# 3P0N True
hTrueTauVisPt3P0N = TH1F('true_3p0n_tau_vis_pt', 'True 3P0N Tau Visible Pt', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt3P0N)

hTrueTauVisTheta3P0N = TH1F('true_3p0n_tau_vis_theta', 'True 3P0N Tau Visible Theta', 10, args.thetaMin, args.thetaMax)
hists.append(hTrueTauVisTheta3P0N)

hTrueTauVisPhi3P0N = TH1F('true_3p0n_tau_vis_phi', 'True 3P0N Tau Visible Phi', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi3P0N)

hTrueTauVisPt3P0N_Barrel = TH1F('true_3p0n_tau_vis_pt_barrel', 'True 3P0N Tau Visible Pt (Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt3P0N_Barrel)

hTrueTauVisPhi3P0N_Barrel = TH1F('true_3p0n_tau_vis_phi_barrel', 'True 3P0N Tau Visible Phi (Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi3P0N_Barrel)

hTrueTauVisPt3P0N_CentBarrel = TH1F('true_3p0n_tau_vis_pt_centbarrel', 'True 3P0N Tau Visible Pt (Central Barrel)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt3P0N_CentBarrel)

hTrueTauVisPhi3P0N_CentBarrel = TH1F('true_3p0n_tau_vis_phi_centbarrel', 'True 3P0N Tau Visible Phi (Central Barrel)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi3P0N_CentBarrel)

hTrueTauVisPt3P0N_Endcap = TH1F('true_3p0n_tau_vis_pt_endcap', 'True 3P0N Tau Visible Pt (Endcap)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt3P0N_Endcap)

hTrueTauVisPhi3P0N_Endcap = TH1F('true_3p0n_tau_vis_phi_endcap', 'True 3P0N Tau Visible Phi (Endcap)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi3P0N_Endcap)

hTrueTauVisPt3P0N_Transition = TH1F('true_3p0n_tau_vis_pt_transition', 'True 3P0N Tau Visible Pt (Transition)', 10, args.pTMin, args.pTMax)
hists.append(hTrueTauVisPt3P0N_Transition)

hTrueTauVisPhi3P0N_Transition = TH1F('true_3p0n_tau_vis_phi_transition', 'True 3P0N Tau Visible Phi (Transition)', 10, args.phiMin, args.phiMax)
hists.append(hTrueTauVisPhi3P0N_Transition)

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

        # Initialize true tau decay mode
        decay_mode = -1

        # Initialize visible observables
        vis_P = []
        vis_pt = 0
        vis_theta = 0
        vis_phi = 0

        # Loop through mc particles
        for mc_particle in mc_particles:

            # Tag taus
            if (abs(mc_particle.getPDG()) == 15):
                
                # Get decay mode
                decay_mode = getDecayMode(mc_particle)

                # Get visible properties
                vis_P = getVisibleP(mc_particle)
                vis_pt = getPt(vis_P)
                vis_theta = getTheta(vis_P)
                vis_phi = getPhi(vis_P)

                # Fill true 1P0N hists
                if (decay_mode == 0):
                    n_1p0n_true += 1
                    hTrueTauVisPt1P0N.Fill(vis_pt)
                    hTrueTauVisTheta1P0N.Fill(vis_theta)
                    hTrueTauVisPhi1P0N.Fill(vis_phi)

                    if (vis_theta > 0.70 and vis_theta < 2.45):
                        hTrueTauVisPt1P0N_Barrel.Fill(vis_pt)
                        hTrueTauVisPhi1P0N_Barrel.Fill(vis_phi)

                    if (vis_theta > 1 and vis_theta < 2):
                        hTrueTauVisPt1P0N_CentBarrel.Fill(vis_pt)
                        hTrueTauVisPhi1P0N_CentBarrel.Fill(vis_phi)

                    elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                        hTrueTauVisPt1P0N_Transition.Fill(vis_pt)
                        hTrueTauVisPhi1P0N_Transition.Fill(vis_phi)

                    elif (vis_theta < 0.577 or vis_theta > 2.56):
                        hTrueTauVisPt1P0N_Endcap.Fill(vis_pt)
                        hTrueTauVisPhi1P0N_Endcap.Fill(vis_phi)

                # Fill true 1P1N hists
                elif (decay_mode == 1):
                    hTrueTauVisPt1P1N.Fill(vis_pt)
                    hTrueTauVisTheta1P1N.Fill(vis_theta)
                    hTrueTauVisPhi1P1N.Fill(vis_phi)

                    if (vis_theta > 0.70 and vis_theta < 2.45):
                        hTrueTauVisPt1P1N_Barrel.Fill(vis_pt)
                        hTrueTauVisPhi1P1N_Barrel.Fill(vis_phi)

                    if (vis_theta > 1 and vis_theta < 2):
                        hTrueTauVisPt1P1N_CentBarrel.Fill(vis_pt)
                        hTrueTauVisPhi1P1N_CentBarrel.Fill(vis_phi)

                    elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                        hTrueTauVisPt1P1N_Transition.Fill(vis_pt)
                        hTrueTauVisPhi1P1N_Transition.Fill(vis_phi)

                    elif (vis_theta < 0.577 or vis_theta > 2.56):
                        hTrueTauVisPt1P1N_Endcap.Fill(vis_pt)
                        hTrueTauVisPhi1P1N_Endcap.Fill(vis_phi)

                # Fill true 1P2N hists
                elif (decay_mode == 2):
                    hTrueTauVisPt1P2N.Fill(vis_pt)
                    hTrueTauVisTheta1P2N.Fill(vis_theta)
                    hTrueTauVisPhi1P2N.Fill(vis_phi)

                    if (vis_theta > 0.70 and vis_theta < 2.45):
                        hTrueTauVisPt1P2N_Barrel.Fill(vis_pt)
                        hTrueTauVisPhi1P2N_Barrel.Fill(vis_phi)

                    if (vis_theta > 1 and vis_theta < 2):
                        hTrueTauVisPt1P2N_CentBarrel.Fill(vis_pt)
                        hTrueTauVisPhi1P2N_CentBarrel.Fill(vis_phi)

                    elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                        hTrueTauVisPt1P2N_Transition.Fill(vis_pt)
                        hTrueTauVisPhi1P2N_Transition.Fill(vis_phi)

                    elif (vis_theta < 0.577 or vis_theta > 2.56):
                        hTrueTauVisPt1P2N_Endcap.Fill(vis_pt)
                        hTrueTauVisPhi1P2N_Endcap.Fill(vis_phi)

                # Fill true 3P0N hists
                elif (decay_mode == 4):
                    n_3p0n_true += 1
                    hTrueTauVisPt3P0N.Fill(vis_pt)
                    hTrueTauVisTheta3P0N.Fill(vis_theta)
                    hTrueTauVisPhi3P0N.Fill(vis_phi)

                    if (vis_theta > 0.70 and vis_theta < 2.45):
                        hTrueTauVisPt3P0N_Barrel.Fill(vis_pt)
                        hTrueTauVisPhi3P0N_Barrel.Fill(vis_phi)

                    if (vis_theta > 1 and vis_theta < 2):
                        hTrueTauVisPt3P0N_CentBarrel.Fill(vis_pt)
                        hTrueTauVisPhi3P0N_CentBarrel.Fill(vis_phi)

                    elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                        hTrueTauVisPt3P0N_Transition.Fill(vis_pt)
                        hTrueTauVisPhi3P0N_Transition.Fill(vis_phi)

                    elif (vis_theta < 0.577 or vis_theta > 2.56):
                        hTrueTauVisPt3P0N_Endcap.Fill(vis_pt)
                        hTrueTauVisPhi3P0N_Endcap.Fill(vis_phi)

                # Fill true 1PXN hists
                if ((decay_mode == 2) or (decay_mode == 1)):
                    hTrueTauVisPt1PXN.Fill(vis_pt)
                    hTrueTauVisTheta1PXN.Fill(vis_theta)
                    hTrueTauVisPhi1PXN.Fill(vis_phi)

                    if (vis_theta > 0.70 and vis_theta < 2.45):
                        hTrueTauVisPt1PXN_Barrel.Fill(vis_pt)
                        hTrueTauVisPhi1PXN_Barrel.Fill(vis_phi)

                    if (vis_theta > 1 and vis_theta < 2):
                        hTrueTauVisPt1PXN_CentBarrel.Fill(vis_pt)
                        hTrueTauVisPhi1PXN_CentBarrel.Fill(vis_phi)

                    elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                        hTrueTauVisPt1PXN_Transition.Fill(vis_pt)
                        hTrueTauVisPhi1PXN_Transition.Fill(vis_phi)

                    elif (vis_theta < 0.577 or vis_theta > 2.56):
                        hTrueTauVisPt1PXN_Endcap.Fill(vis_pt)
                        hTrueTauVisPhi1PXN_Endcap.Fill(vis_phi)
        
        # Loop through reco taus
        for tau in taus:

            # Get number of matched charged pions
            n_matched_pis = getNRecoChargedPis(tau, mc_particles)

            # Get number of neutral pions
            n_pi0s = getNRecoNeutralPis(tau)

            # Fill matched 1P0N hists
            if (n_matched_pis == 1 and decay_mode == 0):
                n_1p0n_reco += 1
                hMatchedTauPt1P0N.Fill(vis_pt)
                hMatchedTauTheta1P0N.Fill(vis_theta)
                hMatchedTauPhi1P0N.Fill(vis_phi)

                if (vis_theta > 0.70 and vis_theta < 2.45):
                    hMatchedTauPt1P0N_Barrel.Fill(vis_pt)
                    hMatchedTauPhi1P0N_Barrel.Fill(vis_phi)

                if (vis_theta > 1 and vis_theta < 2):
                    hMatchedTauPt1P0N_CentBarrel.Fill(vis_pt)
                    hMatchedTauPhi1P0N_CentBarrel.Fill(vis_phi)

                elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                    hMatchedTauPt1P0N_Transition.Fill(vis_pt)
                    hMatchedTauPhi1P0N_Transition.Fill(vis_phi)

                elif (vis_theta < 0.577 or vis_theta > 2.56):
                    hMatchedTauPt1P0N_Endcap.Fill(vis_pt)
                    hMatchedTauPhi1P0N_Endcap.Fill(vis_phi)

            # Fill matched 1P1N hists
            elif (n_matched_pis == 1 and n_pi0s == 1 and decay_mode == 1):
                hMatchedTauPt1P1N.Fill(vis_pt)
                hMatchedTauTheta1P1N.Fill(vis_theta)
                hMatchedTauPhi1P1N.Fill(vis_phi)

                if (vis_theta > 0.70 and vis_theta < 2.45):
                    hMatchedTauPt1P1N_Barrel.Fill(vis_pt)
                    hMatchedTauPhi1P1N_Barrel.Fill(vis_phi)

                if (vis_theta > 1 and vis_theta < 2):
                    hMatchedTauPt1P1N_CentBarrel.Fill(vis_pt)
                    hMatchedTauPhi1P1N_CentBarrel.Fill(vis_phi)

                elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                    hMatchedTauPt1P1N_Transition.Fill(vis_pt)
                    hMatchedTauPhi1P1N_Transition.Fill(vis_phi)

                elif (vis_theta < 0.577 or vis_theta > 2.56):
                    hMatchedTauPt1P1N_Endcap.Fill(vis_pt)
                    hMatchedTauPhi1P1N_Endcap.Fill(vis_phi)

            # Fill matched 1P2N hists
            elif (n_matched_pis == 1 and n_pi0s == 2 and decay_mode == 2):
                hMatchedTauPt1P2N.Fill(vis_pt)
                hMatchedTauTheta1P2N.Fill(vis_theta)
                hMatchedTauPhi1P2N.Fill(vis_phi)

                if (vis_theta > 0.70 and vis_theta < 2.45):
                    hMatchedTauPt1P2N_Barrel.Fill(vis_pt)
                    hMatchedTauPhi1P2N_Barrel.Fill(vis_phi)

                if (vis_theta > 1 and vis_theta < 2):
                    hMatchedTauPt1P2N_CentBarrel.Fill(vis_pt)
                    hMatchedTauPhi1P2N_CentBarrel.Fill(vis_phi)

                elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                    hMatchedTauPt1P2N_Transition.Fill(vis_pt)
                    hMatchedTauPhi1P2N_Transition.Fill(vis_phi)

                elif (vis_theta < 0.577 or vis_theta > 2.56):
                    hMatchedTauPt1P2N_Endcap.Fill(vis_pt)
                    hMatchedTauPhi1P2N_Endcap.Fill(vis_phi)

            # Fill matched 3P0N hists
            elif (n_matched_pis == 3 and decay_mode == 4):
                n_3p0n_reco += 1
                hMatchedTauPt3P0N.Fill(vis_pt)
                hMatchedTauTheta3P0N.Fill(vis_theta)
                hMatchedTauPhi3P0N.Fill(vis_phi)
                
                if (vis_theta > 0.70 and vis_theta < 2.45):
                    hMatchedTauPt3P0N_Barrel.Fill(vis_pt)
                    hMatchedTauPhi3P0N_Barrel.Fill(vis_phi)

                if (vis_theta > 1 and vis_theta < 2):
                    hMatchedTauPt3P0N_CentBarrel.Fill(vis_pt)
                    hMatchedTauPhi3P0N_CentBarrel.Fill(vis_phi)

                elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                    hMatchedTauPt3P0N_Transition.Fill(vis_pt)
                    hMatchedTauPhi3P0N_Transition.Fill(vis_phi)

                elif (vis_theta < 0.577 or vis_theta > 2.56):
                    hMatchedTauPt3P0N_Endcap.Fill(vis_pt)
                    hMatchedTauPhi3P0N_Endcap.Fill(vis_phi)

            # Fill matched 1PXN hists
            if (n_matched_pis == 1 and n_pi0s > 0 and (decay_mode == 2 or decay_mode == 1)):
                hMatchedTauPt1PXN.Fill(vis_pt)
                hMatchedTauTheta1PXN.Fill(vis_theta)
                hMatchedTauPhi1PXN.Fill(vis_phi)

                if (vis_theta > 0.70 and vis_theta < 2.45):
                    hMatchedTauPt1PXN_Barrel.Fill(vis_pt)
                    hMatchedTauPhi1PXN_Barrel.Fill(vis_phi)

                if (vis_theta > 1 and vis_theta < 2):
                    hMatchedTauPt1PXN_CentBarrel.Fill(vis_pt)
                    hMatchedTauPhi1PXN_CentBarrel.Fill(vis_phi)

                elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                    hMatchedTauPt1PXN_Transition.Fill(vis_pt)
                    hMatchedTauPhi1PXN_Transition.Fill(vis_phi)

                elif (vis_theta < 0.577 or vis_theta > 2.56):
                    hMatchedTauPt1PXN_Endcap.Fill(vis_pt)
                    hMatchedTauPhi1PXN_Endcap.Fill(vis_phi)

    reader.close()

print(f'1p0n eff: {n_1p0n_reco/n_1p0n_true}')
print(f'3p0n eff: {n_3p0n_reco/n_3p0n_true}')

# Create 1P0N efficiency hists
hPtEff1P0N = hMatchedTauPt1P0N.Clone('1p0n_tau_pt_eff')
hPtEff1P0N.Divide(hPtEff1P0N, hTrueTauVisPt1P0N, 1, 1, 'B')
hPtEff1P0N.SetLineColor(6)
hPtEff1P0N.SetLineWidth(2)
hPtEff1P0N.SetTitle('1P0N Reconstruction Efficiency vs Pt')
hPtEff1P0N.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N.SetStats(0)
hists.append(hPtEff1P0N)

hThetaEff1P0N = hMatchedTauTheta1P0N.Clone('1p0n_tau_theta_eff')
hThetaEff1P0N.Divide(hThetaEff1P0N, hTrueTauVisTheta1P0N, 1, 1, 'B')
hThetaEff1P0N.SetLineColor(7)
hThetaEff1P0N.SetLineWidth(2)
hThetaEff1P0N.SetTitle('1P0N Reconstruction Efficiency vs Theta')
hThetaEff1P0N.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hThetaEff1P0N.GetYaxis().SetTitle('#epsilon')
hThetaEff1P0N.SetStats(0)
hists.append(hThetaEff1P0N)

hPhiEff1P0N = hMatchedTauPhi1P0N.Clone('1p0n_tau_phi_eff')
hPhiEff1P0N.Divide(hPhiEff1P0N, hTrueTauVisPhi1P0N, 1, 1, 'B')
hPhiEff1P0N.SetLineColor(418)
hPhiEff1P0N.SetLineWidth(2)
hPhiEff1P0N.SetTitle('1P0N Reconstruction Efficiency vs Phi')
hPhiEff1P0N.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N.SetStats(0)
hists.append(hPhiEff1P0N)

hPtEff1P0N_Barrel = hMatchedTauPt1P0N_Barrel.Clone('1p0n_tau_pt_eff_barrel')
hPtEff1P0N_Barrel.Divide(hPtEff1P0N_Barrel, hTrueTauVisPt1P0N_Barrel, 1, 1, 'B')
hPtEff1P0N_Barrel.SetLineColor(6)
hPtEff1P0N_Barrel.SetLineWidth(2)
hPtEff1P0N_Barrel.SetTitle('1P0N Reconstruction Efficiency vs Pt (Barrel)')
hPtEff1P0N_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N_Barrel.SetStats(0)
hists.append(hPtEff1P0N_Barrel)

hPhiEff1P0N_Barrel = hMatchedTauPhi1P0N_Barrel.Clone('1p0n_tau_phi_eff_barrel')
hPhiEff1P0N_Barrel.Divide(hPhiEff1P0N_Barrel, hTrueTauVisPhi1P0N_Barrel, 1, 1, 'B')
hPhiEff1P0N_Barrel.SetLineColor(418)
hPhiEff1P0N_Barrel.SetLineWidth(2)
hPhiEff1P0N_Barrel.SetTitle('1P0N Reconstruction Efficiency vs Phi (Barrel)')
hPhiEff1P0N_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N_Barrel.SetStats(0)
hists.append(hPhiEff1P0N_Barrel)

hPtEff1P0N_CentBarrel = hMatchedTauPt1P0N_CentBarrel.Clone('1p0n_tau_pt_eff_centbarrel')
hPtEff1P0N_CentBarrel.Divide(hPtEff1P0N_CentBarrel, hTrueTauVisPt1P0N_CentBarrel, 1, 1, 'B')
hPtEff1P0N_CentBarrel.SetLineColor(6)
hPtEff1P0N_CentBarrel.SetLineWidth(2)
hPtEff1P0N_CentBarrel.SetTitle('1P0N Reconstruction Efficiency vs Pt (Central Barrel)')
hPtEff1P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N_CentBarrel.SetStats(0)
hists.append(hPtEff1P0N_CentBarrel)

hPhiEff1P0N_CentBarrel = hMatchedTauPhi1P0N_CentBarrel.Clone('1p0n_tau_phi_eff_centbarrel')
hPhiEff1P0N_CentBarrel.Divide(hPhiEff1P0N_CentBarrel, hTrueTauVisPhi1P0N_CentBarrel, 1, 1, 'B')
hPhiEff1P0N_CentBarrel.SetLineColor(418)
hPhiEff1P0N_CentBarrel.SetLineWidth(2)
hPhiEff1P0N_CentBarrel.SetTitle('1P0N Reconstruction Efficiency vs Phi (Central Barrel)')
hPhiEff1P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N_CentBarrel.SetStats(0)
hists.append(hPhiEff1P0N_CentBarrel)

hPtEff1P0N_Endcap = hMatchedTauPt1P0N_Endcap.Clone('1p0n_tau_pt_eff_endcap')
hPtEff1P0N_Endcap.Divide(hPtEff1P0N_Endcap, hTrueTauVisPt1P0N_Endcap, 1, 1, 'B')
hPtEff1P0N_Endcap.SetLineColor(6)
hPtEff1P0N_Endcap.SetLineWidth(2)
hPtEff1P0N_Endcap.SetTitle('1P0N Reconstruction Efficiency vs Pt (Endcap)')
hPtEff1P0N_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N_Endcap.SetStats(0)
hists.append(hPtEff1P0N_Endcap)

hPhiEff1P0N_Endcap = hMatchedTauPhi1P0N_Endcap.Clone('1p0n_tau_phi_eff_endcap')
hPhiEff1P0N_Endcap.Divide(hPhiEff1P0N_Endcap, hTrueTauVisPhi1P0N_Endcap, 1, 1, 'B')
hPhiEff1P0N_Endcap.SetLineColor(418)
hPhiEff1P0N_Endcap.SetLineWidth(2)
hPhiEff1P0N_Endcap.SetTitle('1P0N Reconstruction Efficiency vs Phi (Endcap)')
hPhiEff1P0N_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N_Endcap.SetStats(0)
hists.append(hPhiEff1P0N_Endcap)

hPtEff1P0N_Transition = hMatchedTauPt1P0N_Transition.Clone('1p0n_tau_pt_eff_transition')
hPtEff1P0N_Transition.Divide(hPtEff1P0N_Transition, hTrueTauVisPt1P0N_Transition, 1, 1, 'B')
hPtEff1P0N_Transition.SetLineColor(6)
hPtEff1P0N_Transition.SetLineWidth(2)
hPtEff1P0N_Transition.SetTitle('1P0N Reconstruction Efficiency vs Pt (Transition)')
hPtEff1P0N_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N_Transition.SetStats(0)
hists.append(hPtEff1P0N_Transition)

hPhiEff1P0N_Transition = hMatchedTauPhi1P0N_Transition.Clone('1p0n_tau_phi_eff_transition')
hPhiEff1P0N_Transition.Divide(hPhiEff1P0N_Transition, hTrueTauVisPhi1P0N_Transition, 1, 1, 'B')
hPhiEff1P0N_Transition.SetLineColor(418)
hPhiEff1P0N_Transition.SetLineWidth(2)
hPhiEff1P0N_Transition.SetTitle('1P0N Reconstruction Efficiency vs Phi (Transition)')
hPhiEff1P0N_Transition.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N_Transition.SetStats(0)
hists.append(hPhiEff1P0N_Transition)

# Create 1P1N efficiency hists
hPtEff1P1N = hMatchedTauPt1P1N.Clone('1p1n_tau_pt_eff')
hPtEff1P1N.Divide(hPtEff1P1N, hTrueTauVisPt1P1N, 1, 1, 'B')
hPtEff1P1N.SetLineColor(6)
hPtEff1P1N.SetLineWidth(2)
hPtEff1P1N.SetTitle('1P1N Reconstruction Efficiency vs Pt')
hPtEff1P1N.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P1N.GetYaxis().SetTitle('#epsilon')
hPtEff1P1N.SetStats(1)
hists.append(hPtEff1P1N)

hThetaEff1P1N = hMatchedTauTheta1P1N.Clone('1p1n_tau_theta_eff')
hThetaEff1P1N.Divide(hThetaEff1P1N, hTrueTauVisTheta1P1N, 1, 1, 'B')
hThetaEff1P1N.SetLineColor(7)
hThetaEff1P1N.SetLineWidth(2)
hThetaEff1P1N.SetTitle('1P1N Reconstruction Efficiency vs Theta')
hThetaEff1P1N.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hThetaEff1P1N.GetYaxis().SetTitle('#epsilon')
hThetaEff1P1N.SetStats(1)
hists.append(hThetaEff1P1N)

hPhiEff1P1N = hMatchedTauPhi1P1N.Clone('1p1n_tau_phi_eff')
hPhiEff1P1N.Divide(hPhiEff1P1N, hTrueTauVisPhi1P1N, 1, 1, 'B')
hPhiEff1P1N.SetLineColor(418)
hPhiEff1P1N.SetLineWidth(2)
hPhiEff1P1N.SetTitle('1P1N Reconstruction Efficiency vs Phi')
hPhiEff1P1N.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P1N.GetYaxis().SetTitle('#epsilon')
hPhiEff1P1N.SetStats(1)
hists.append(hPhiEff1P1N)

hPtEff1P1N_Barrel = hMatchedTauPt1P1N_Barrel.Clone('1p1n_tau_pt_eff_barrel')
hPtEff1P1N_Barrel.Divide(hPtEff1P1N_Barrel, hTrueTauVisPt1P1N_Barrel, 1, 1, 'B')
hPtEff1P1N_Barrel.SetLineColor(6)
hPtEff1P1N_Barrel.SetLineWidth(2)
hPtEff1P1N_Barrel.SetTitle('1P1N Reconstruction Efficiency vs Pt (Barrel)')
hPtEff1P1N_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P1N_Barrel.GetYaxis().SetTitle('#epsilon')
hPtEff1P1N_Barrel.SetStats(1)
hists.append(hPtEff1P1N_Barrel)

hPhiEff1P1N_Barrel = hMatchedTauPhi1P1N_Barrel.Clone('1p1n_tau_phi_eff_barrel')
hPhiEff1P1N_Barrel.Divide(hPhiEff1P1N_Barrel, hTrueTauVisPhi1P1N_Barrel, 1, 1, 'B')
hPhiEff1P1N_Barrel.SetLineColor(418)
hPhiEff1P1N_Barrel.SetLineWidth(2)
hPhiEff1P1N_Barrel.SetTitle('1P1N Reconstruction Efficiency vs Phi (Barrel)')
hPhiEff1P1N_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P1N_Barrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1P1N_Barrel.SetStats(1)
hists.append(hPhiEff1P1N_Barrel)

hPtEff1P1N_CentBarrel = hMatchedTauPt1P1N_CentBarrel.Clone('1p1n_tau_pt_eff_centbarrel')
hPtEff1P1N_CentBarrel.Divide(hPtEff1P1N_CentBarrel, hTrueTauVisPt1P1N_CentBarrel, 1, 1, 'B')
hPtEff1P1N_CentBarrel.SetLineColor(6)
hPtEff1P1N_CentBarrel.SetLineWidth(2)
hPtEff1P1N_CentBarrel.SetTitle('1P1N Reconstruction Efficiency vs Pt (Central Barrel)')
hPtEff1P1N_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P1N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPtEff1P1N_CentBarrel.SetStats(1)
hists.append(hPtEff1P1N_CentBarrel)

hPhiEff1P1N_CentBarrel = hMatchedTauPhi1P1N_CentBarrel.Clone('1p1n_tau_phi_eff_centbarrel')
hPhiEff1P1N_CentBarrel.Divide(hPhiEff1P1N_CentBarrel, hTrueTauVisPhi1P1N_CentBarrel, 1, 1, 'B')
hPhiEff1P1N_CentBarrel.SetLineColor(418)
hPhiEff1P1N_CentBarrel.SetLineWidth(2)
hPhiEff1P1N_CentBarrel.SetTitle('1P1N Reconstruction Efficiency vs Phi (Central Barrel)')
hPhiEff1P1N_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P1N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1P1N_CentBarrel.SetStats(1)
hists.append(hPhiEff1P1N_CentBarrel)

hPtEff1P1N_Endcap = hMatchedTauPt1P1N_Endcap.Clone('1p1n_tau_pt_eff_endcap')
hPtEff1P1N_Endcap.Divide(hPtEff1P1N_Endcap, hTrueTauVisPt1P1N_Endcap, 1, 1, 'B')
hPtEff1P1N_Endcap.SetLineColor(6)
hPtEff1P1N_Endcap.SetLineWidth(2)
hPtEff1P1N_Endcap.SetTitle('1P1N Reconstruction Efficiency vs Pt (Endcap)')
hPtEff1P1N_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P1N_Endcap.GetYaxis().SetTitle('#epsilon')
hPtEff1P1N_Endcap.SetStats(1)
hists.append(hPtEff1P1N_Endcap)

hPhiEff1P1N_Endcap = hMatchedTauPhi1P1N_Endcap.Clone('1p1n_tau_phi_eff_endcap')
hPhiEff1P1N_Endcap.Divide(hPhiEff1P1N_Endcap, hTrueTauVisPhi1P1N_Endcap, 1, 1, 'B')
hPhiEff1P1N_Endcap.SetLineColor(418)
hPhiEff1P1N_Endcap.SetLineWidth(2)
hPhiEff1P1N_Endcap.SetTitle('1P1N Reconstruction Efficiency vs Phi (Endcap)')
hPhiEff1P1N_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P1N_Endcap.GetYaxis().SetTitle('#epsilon')
hPhiEff1P1N_Endcap.SetStats(1)
hists.append(hPhiEff1P1N_Endcap)

hPtEff1P1N_Transition = hMatchedTauPt1P1N_Transition.Clone('1p1n_tau_pt_eff_transition')
hPtEff1P1N_Transition.Divide(hPtEff1P1N_Transition, hTrueTauVisPt1P1N_Transition, 1, 1, 'B')
hPtEff1P1N_Transition.SetLineColor(6)
hPtEff1P1N_Transition.SetLineWidth(2)
hPtEff1P1N_Transition.SetTitle('1P1N Reconstruction Efficiency vs Pt (Transition)')
hPtEff1P1N_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P1N_Transition.GetYaxis().SetTitle('#epsilon')
hPtEff1P1N_Transition.SetStats(1)
hists.append(hPtEff1P1N_Transition)

hPhiEff1P1N_Transition = hMatchedTauPhi1P1N_Transition.Clone('1p1n_tau_phi_eff_transition')
hPhiEff1P1N_Transition.Divide(hPhiEff1P1N_Transition, hTrueTauVisPhi1P1N_Transition, 1, 1, 'B')
hPhiEff1P1N_Transition.SetLineColor(418)
hPhiEff1P1N_Transition.SetLineWidth(2)
hPhiEff1P1N_Transition.SetTitle('1P1N Reconstruction Efficiency vs Phi (Transition)')
hPhiEff1P1N_Transition.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P1N_Transition.GetYaxis().SetTitle('#epsilon')
hPhiEff1P1N_Transition.SetStats(1)
hists.append(hPhiEff1P1N_Transition)

# Create 1P2N efficiency hists
hPtEff1P2N = hMatchedTauPt1P2N.Clone('1p2n_tau_pt_eff')
hPtEff1P2N.Divide(hPtEff1P2N, hTrueTauVisPt1P2N, 1, 1, 'B')
hPtEff1P2N.SetLineColor(6)
hPtEff1P2N.SetLineWidth(2)
hPtEff1P2N.SetTitle('1P2N Reconstruction Efficiency vs Pt')
hPtEff1P2N.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P2N.GetYaxis().SetTitle('#epsilon')
hPtEff1P2N.SetStats(1)
hists.append(hPtEff1P2N)

hThetaEff1P2N = hMatchedTauTheta1P2N.Clone('1p2n_tau_theta_eff')
hThetaEff1P2N.Divide(hThetaEff1P2N, hTrueTauVisTheta1P2N, 1, 1, 'B')
hThetaEff1P2N.SetLineColor(7)
hThetaEff1P2N.SetLineWidth(2)
hThetaEff1P2N.SetTitle('1P2N Reconstruction Efficiency vs Theta')
hThetaEff1P2N.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hThetaEff1P2N.GetYaxis().SetTitle('#epsilon')
hThetaEff1P2N.SetStats(1)
hists.append(hThetaEff1P2N)

hPhiEff1P2N = hMatchedTauPhi1P2N.Clone('1p2n_tau_phi_eff')
hPhiEff1P2N.Divide(hPhiEff1P2N, hTrueTauVisPhi1P2N, 1, 1, 'B')
hPhiEff1P2N.SetLineColor(418)
hPhiEff1P2N.SetLineWidth(2)
hPhiEff1P2N.SetTitle('1P2N Reconstruction Efficiency vs Phi')
hPhiEff1P2N.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P2N.GetYaxis().SetTitle('#epsilon')
hPhiEff1P2N.SetStats(1)
hists.append(hPhiEff1P2N)

hPtEff1P2N_Barrel = hMatchedTauPt1P2N_Barrel.Clone('1p2n_tau_pt_eff_barrel')
hPtEff1P2N_Barrel.Divide(hPtEff1P2N_Barrel, hTrueTauVisPt1P2N_Barrel, 1, 1, 'B')
hPtEff1P2N_Barrel.SetLineColor(6)
hPtEff1P2N_Barrel.SetLineWidth(2)
hPtEff1P2N_Barrel.SetTitle('1P2N Reconstruction Efficiency vs Pt (Barrel)')
hPtEff1P2N_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P2N_Barrel.GetYaxis().SetTitle('#epsilon')
hPtEff1P2N_Barrel.SetStats(1)
hists.append(hPtEff1P2N_Barrel)

hPhiEff1P2N_Barrel = hMatchedTauPhi1P2N_Barrel.Clone('1p2n_tau_phi_eff_barrel')
hPhiEff1P2N_Barrel.Divide(hPhiEff1P2N_Barrel, hTrueTauVisPhi1P2N_Barrel, 1, 1, 'B')
hPhiEff1P2N_Barrel.SetLineColor(418)
hPhiEff1P2N_Barrel.SetLineWidth(2)
hPhiEff1P2N_Barrel.SetTitle('1P2N Reconstruction Efficiency vs Phi (Barrel)')
hPhiEff1P2N_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P2N_Barrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1P2N_Barrel.SetStats(1)
hists.append(hPhiEff1P2N_Barrel)

hPtEff1P2N_CentBarrel = hMatchedTauPt1P2N_CentBarrel.Clone('1p2n_tau_pt_eff_centbarrel')
hPtEff1P2N_CentBarrel.Divide(hPtEff1P2N_CentBarrel, hTrueTauVisPt1P2N_CentBarrel, 1, 1, 'B')
hPtEff1P2N_CentBarrel.SetLineColor(6)
hPtEff1P2N_CentBarrel.SetLineWidth(2)
hPtEff1P2N_CentBarrel.SetTitle('1P2N Reconstruction Efficiency vs Pt (Central Barrel)')
hPtEff1P2N_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P2N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPtEff1P2N_CentBarrel.SetStats(1)
hists.append(hPtEff1P2N_CentBarrel)

hPhiEff1P2N_CentBarrel = hMatchedTauPhi1P2N_CentBarrel.Clone('1p2n_tau_phi_eff_centbarrel')
hPhiEff1P2N_CentBarrel.Divide(hPhiEff1P2N_CentBarrel, hTrueTauVisPhi1P2N_CentBarrel, 1, 1, 'B')
hPhiEff1P2N_CentBarrel.SetLineColor(418)
hPhiEff1P2N_CentBarrel.SetLineWidth(2)
hPhiEff1P2N_CentBarrel.SetTitle('1P2N Reconstruction Efficiency vs Phi (Central Barrel)')
hPhiEff1P2N_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P2N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1P2N_CentBarrel.SetStats(1)
hists.append(hPhiEff1P2N_CentBarrel)

hPtEff1P2N_Endcap = hMatchedTauPt1P2N_Endcap.Clone('1p2n_tau_pt_eff_endcap')
hPtEff1P2N_Endcap.Divide(hPtEff1P2N_Endcap, hTrueTauVisPt1P2N_Endcap, 1, 1, 'B')
hPtEff1P2N_Endcap.SetLineColor(6)
hPtEff1P2N_Endcap.SetLineWidth(2)
hPtEff1P2N_Endcap.SetTitle('1P2N Reconstruction Efficiency vs Pt (Endcap)')
hPtEff1P2N_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P2N_Endcap.GetYaxis().SetTitle('#epsilon')
hPtEff1P2N_Endcap.SetStats(1)
hists.append(hPtEff1P2N_Endcap)

hPhiEff1P2N_Endcap = hMatchedTauPhi1P2N_Endcap.Clone('1p2n_tau_phi_eff_endcap')
hPhiEff1P2N_Endcap.Divide(hPhiEff1P2N_Endcap, hTrueTauVisPhi1P2N_Endcap, 1, 1, 'B')
hPhiEff1P2N_Endcap.SetLineColor(418)
hPhiEff1P2N_Endcap.SetLineWidth(2)
hPhiEff1P2N_Endcap.SetTitle('1P2N Reconstruction Efficiency vs Phi (Endcap)')
hPhiEff1P2N_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P2N_Endcap.GetYaxis().SetTitle('#epsilon')
hPhiEff1P2N_Endcap.SetStats(1)
hists.append(hPhiEff1P2N_Endcap)

hPtEff1P2N_Transition = hMatchedTauPt1P2N_Transition.Clone('1p2n_tau_pt_eff_transition')
hPtEff1P2N_Transition.Divide(hPtEff1P2N_Transition, hTrueTauVisPt1P2N_Transition, 1, 1, 'B')
hPtEff1P2N_Transition.SetLineColor(6)
hPtEff1P2N_Transition.SetLineWidth(2)
hPtEff1P2N_Transition.SetTitle('1P2N Reconstruction Efficiency vs Pt (Transition)')
hPtEff1P2N_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P2N_Transition.GetYaxis().SetTitle('#epsilon')
hPtEff1P2N_Transition.SetStats(1)
hists.append(hPtEff1P2N_Transition)

hPhiEff1P2N_Transition = hMatchedTauPhi1P2N_Transition.Clone('1p2n_tau_phi_eff_transition')
hPhiEff1P2N_Transition.Divide(hPhiEff1P2N_Transition, hTrueTauVisPhi1P2N_Transition, 1, 1, 'B')
hPhiEff1P2N_Transition.SetLineColor(418)
hPhiEff1P2N_Transition.SetLineWidth(2)
hPhiEff1P2N_Transition.SetTitle('1P2N Reconstruction Efficiency vs Phi (Transition)')
hPhiEff1P2N_Transition.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P2N_Transition.GetYaxis().SetTitle('#epsilon')
hPhiEff1P2N_Transition.SetStats(1)
hists.append(hPhiEff1P2N_Transition)

# Create 1PXN efficiency hists
hPtEff1PXN = hMatchedTauPt1PXN.Clone('1pXn_tau_pt_eff')
hPtEff1PXN.Divide(hPtEff1PXN, hTrueTauVisPt1PXN, 1, 1, 'B')
hPtEff1PXN.SetLineColor(6)
hPtEff1PXN.SetLineWidth(2)
hPtEff1PXN.SetTitle('1PXN Reconstruction Efficiency vs Pt')
hPtEff1PXN.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1PXN.GetYaxis().SetTitle('#epsilon')
hPtEff1PXN.SetStats(1)
hists.append(hPtEff1PXN)

hThetaEff1PXN = hMatchedTauTheta1PXN.Clone('1pXn_tau_theta_eff')
hThetaEff1PXN.Divide(hThetaEff1PXN, hTrueTauVisTheta1PXN, 1, 1, 'B')
hThetaEff1PXN.SetLineColor(7)
hThetaEff1PXN.SetLineWidth(2)
hThetaEff1PXN.SetTitle('1PXN Reconstruction Efficiency vs Theta')
hThetaEff1PXN.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hThetaEff1PXN.GetYaxis().SetTitle('#epsilon')
hThetaEff1PXN.SetStats(1)
hists.append(hThetaEff1PXN)

hPhiEff1PXN = hMatchedTauPhi1PXN.Clone('1pXn_tau_phi_eff')
hPhiEff1PXN.Divide(hPhiEff1PXN, hTrueTauVisPhi1PXN, 1, 1, 'B')
hPhiEff1PXN.SetLineColor(418)
hPhiEff1PXN.SetLineWidth(2)
hPhiEff1PXN.SetTitle('1PXN Reconstruction Efficiency vs Phi')
hPhiEff1PXN.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1PXN.GetYaxis().SetTitle('#epsilon')
hPhiEff1PXN.SetStats(1)
hists.append(hPhiEff1PXN)

hPtEff1PXN_Barrel = hMatchedTauPt1PXN_Barrel.Clone('1pXn_tau_pt_eff_barrel')
hPtEff1PXN_Barrel.Divide(hPtEff1PXN_Barrel, hTrueTauVisPt1PXN_Barrel, 1, 1, 'B')
hPtEff1PXN_Barrel.SetLineColor(6)
hPtEff1PXN_Barrel.SetLineWidth(2)
hPtEff1PXN_Barrel.SetTitle('1PXN Reconstruction Efficiency vs Pt (Barrel)')
hPtEff1PXN_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1PXN_Barrel.GetYaxis().SetTitle('#epsilon')
hPtEff1PXN_Barrel.SetStats(1)
hists.append(hPtEff1PXN_Barrel)

hPhiEff1PXN_Barrel = hMatchedTauPhi1PXN_Barrel.Clone('1pXn_tau_phi_eff_barrel')
hPhiEff1PXN_Barrel.Divide(hPhiEff1PXN_Barrel, hTrueTauVisPhi1PXN_Barrel, 1, 1, 'B')
hPhiEff1PXN_Barrel.SetLineColor(418)
hPhiEff1PXN_Barrel.SetLineWidth(2)
hPhiEff1PXN_Barrel.SetTitle('1PXN Reconstruction Efficiency vs Phi (Barrel)')
hPhiEff1PXN_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1PXN_Barrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1PXN_Barrel.SetStats(1)
hists.append(hPhiEff1PXN_Barrel)

hPtEff1PXN_CentBarrel = hMatchedTauPt1PXN_CentBarrel.Clone('1pXn_tau_pt_eff_centbarrel')
hPtEff1PXN_CentBarrel.Divide(hPtEff1PXN_CentBarrel, hTrueTauVisPt1PXN_CentBarrel, 1, 1, 'B')
hPtEff1PXN_CentBarrel.SetLineColor(6)
hPtEff1PXN_CentBarrel.SetLineWidth(2)
hPtEff1PXN_CentBarrel.SetTitle('1PXN Reconstruction Efficiency vs Pt (Central Barrel)')
hPtEff1PXN_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1PXN_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPtEff1PXN_CentBarrel.SetStats(1)
hists.append(hPtEff1PXN_CentBarrel)

hPhiEff1PXN_CentBarrel = hMatchedTauPhi1PXN_CentBarrel.Clone('1pXn_tau_phi_eff_centbarrel')
hPhiEff1PXN_CentBarrel.Divide(hPhiEff1PXN_CentBarrel, hTrueTauVisPhi1PXN_CentBarrel, 1, 1, 'B')
hPhiEff1PXN_CentBarrel.SetLineColor(418)
hPhiEff1PXN_CentBarrel.SetLineWidth(2)
hPhiEff1PXN_CentBarrel.SetTitle('1PXN Reconstruction Efficiency vs Phi (Central Barrel)')
hPhiEff1PXN_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1PXN_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1PXN_CentBarrel.SetStats(1)
hists.append(hPhiEff1PXN_CentBarrel)

hPtEff1PXN_Endcap = hMatchedTauPt1PXN_Endcap.Clone('1pXn_tau_pt_eff_endcap')
hPtEff1PXN_Endcap.Divide(hPtEff1PXN_Endcap, hTrueTauVisPt1PXN_Endcap, 1, 1, 'B')
hPtEff1PXN_Endcap.SetLineColor(6)
hPtEff1PXN_Endcap.SetLineWidth(2)
hPtEff1PXN_Endcap.SetTitle('1PXN Reconstruction Efficiency vs Pt (Endcap)')
hPtEff1PXN_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1PXN_Endcap.GetYaxis().SetTitle('#epsilon')
hPtEff1PXN_Endcap.SetStats(1)
hists.append(hPtEff1PXN_Endcap)

hPhiEff1PXN_Endcap = hMatchedTauPhi1PXN_Endcap.Clone('1pXn_tau_phi_eff_endcap')
hPhiEff1PXN_Endcap.Divide(hPhiEff1PXN_Endcap, hTrueTauVisPhi1PXN_Endcap, 1, 1, 'B')
hPhiEff1PXN_Endcap.SetLineColor(418)
hPhiEff1PXN_Endcap.SetLineWidth(2)
hPhiEff1PXN_Endcap.SetTitle('1PXN Reconstruction Efficiency vs Phi (Endcap)')
hPhiEff1PXN_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1PXN_Endcap.GetYaxis().SetTitle('#epsilon')
hPhiEff1PXN_Endcap.SetStats(1)
hists.append(hPhiEff1PXN_Endcap)

hPtEff1PXN_Transition = hMatchedTauPt1PXN_Transition.Clone('1pXn_tau_pt_eff_transition')
hPtEff1PXN_Transition.Divide(hPtEff1PXN_Transition, hTrueTauVisPt1PXN_Transition, 1, 1, 'B')
hPtEff1PXN_Transition.SetLineColor(6)
hPtEff1PXN_Transition.SetLineWidth(2)
hPtEff1PXN_Transition.SetTitle('1PXN Reconstruction Efficiency vs Pt (Transition)')
hPtEff1PXN_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1PXN_Transition.GetYaxis().SetTitle('#epsilon')
hPtEff1PXN_Transition.SetStats(1)
hists.append(hPtEff1PXN_Transition)

hPhiEff1PXN_Transition = hMatchedTauPhi1PXN_Transition.Clone('1pXn_tau_phi_eff_transition')
hPhiEff1PXN_Transition.Divide(hPhiEff1PXN_Transition, hTrueTauVisPhi1PXN_Transition, 1, 1, 'B')
hPhiEff1PXN_Transition.SetLineColor(418)
hPhiEff1PXN_Transition.SetLineWidth(2)
hPhiEff1PXN_Transition.SetTitle('1PXN Reconstruction Efficiency vs Phi (Transition)')
hPhiEff1PXN_Transition.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1PXN_Transition.GetYaxis().SetTitle('#epsilon')
hPhiEff1PXN_Transition.SetStats(1)
hists.append(hPhiEff1PXN_Transition)

# Create 3P0N efficiency hists
hPtEff3P0N = hMatchedTauPt3P0N.Clone('3p0n_tau_pt_eff')
hPtEff3P0N.Divide(hPtEff3P0N, hTrueTauVisPt3P0N, 1, 1, 'B')
hPtEff3P0N.SetLineColor(6)
hPtEff3P0N.SetLineWidth(2)
hPtEff3P0N.SetTitle('3P0N Reconstruction Efficiency vs Pt')
hPtEff3P0N.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N.SetStats(0)
hists.append(hPtEff3P0N)

hThetaEff3P0N = hMatchedTauTheta3P0N.Clone('3p0n_tau_theta_eff')
hThetaEff3P0N.Divide(hThetaEff3P0N, hTrueTauVisTheta3P0N, 1, 1, 'B')
hThetaEff3P0N.SetLineColor(7)
hThetaEff3P0N.SetLineWidth(2)
hThetaEff3P0N.SetTitle('3P0N Reconstruction Efficiency vs Theta')
hThetaEff3P0N.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hThetaEff3P0N.GetYaxis().SetTitle('#epsilon')
hThetaEff3P0N.SetStats(0)
hists.append(hThetaEff3P0N)

hPhiEff3P0N = hMatchedTauPhi3P0N.Clone('3p0n_tau_phi_eff')
hPhiEff3P0N.Divide(hPhiEff3P0N, hTrueTauVisPhi3P0N, 1, 1, 'B')
hPhiEff3P0N.SetLineColor(418)
hPhiEff3P0N.SetLineWidth(2)
hPhiEff3P0N.SetTitle('3P0N Reconstruction Efficiency vs Phi')
hPhiEff3P0N.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N.SetStats(0)
hists.append(hPhiEff3P0N)

hPtEff3P0N_Barrel = hMatchedTauPt3P0N_Barrel.Clone('3p0n_tau_pt_eff_barrel')
hPtEff3P0N_Barrel.Divide(hPtEff3P0N_Barrel, hTrueTauVisPt3P0N_Barrel, 1, 1, 'B')
hPtEff3P0N_Barrel.SetLineColor(6)
hPtEff3P0N_Barrel.SetLineWidth(2)
hPtEff3P0N_Barrel.SetTitle('3P0N Reconstruction Efficiency vs Pt (Barrel)')
hPtEff3P0N_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N_Barrel.SetStats(0)
hists.append(hPtEff3P0N_Barrel)

hPhiEff3P0N_Barrel = hMatchedTauPhi3P0N_Barrel.Clone('3p0n_tau_phi_eff_barrel')
hPhiEff3P0N_Barrel.Divide(hPhiEff3P0N_Barrel, hTrueTauVisPhi3P0N_Barrel, 1, 1, 'B')
hPhiEff3P0N_Barrel.SetLineColor(418)
hPhiEff3P0N_Barrel.SetLineWidth(2)
hPhiEff3P0N_Barrel.SetTitle('3P0N Reconstruction Efficiency vs Phi (Barrel)')
hPhiEff3P0N_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N_Barrel.SetStats(0)
hists.append(hPhiEff3P0N_Barrel)

hPtEff3P0N_CentBarrel = hMatchedTauPt3P0N_CentBarrel.Clone('3p0n_tau_pt_eff_centbarrel')
hPtEff3P0N_CentBarrel.Divide(hPtEff3P0N_CentBarrel, hTrueTauVisPt3P0N_CentBarrel, 1, 1, 'B')
hPtEff3P0N_CentBarrel.SetLineColor(6)
hPtEff3P0N_CentBarrel.SetLineWidth(2)
hPtEff3P0N_CentBarrel.SetTitle('3P0N Reconstruction Efficiency vs Pt (Central Barrel)')
hPtEff3P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N_CentBarrel.SetStats(0)
hists.append(hPtEff3P0N_CentBarrel)

hPhiEff3P0N_CentBarrel = hMatchedTauPhi3P0N_CentBarrel.Clone('3p0n_tau_phi_eff_centbarrel')
hPhiEff3P0N_CentBarrel.Divide(hPhiEff3P0N_CentBarrel, hTrueTauVisPhi3P0N_CentBarrel, 1, 1, 'B')
hPhiEff3P0N_CentBarrel.SetLineColor(418)
hPhiEff3P0N_CentBarrel.SetLineWidth(2)
hPhiEff3P0N_CentBarrel.SetTitle('3P0N Reconstruction Efficiency vs Phi (Central Barrel)')
hPhiEff3P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N_CentBarrel.SetStats(0)
hists.append(hPhiEff3P0N_CentBarrel)

hPtEff3P0N_Endcap = hMatchedTauPt3P0N_Endcap.Clone('3p0n_tau_pt_eff_endcap')
hPtEff3P0N_Endcap.Divide(hPtEff3P0N_Endcap, hTrueTauVisPt3P0N_Endcap, 1, 1, 'B')
hPtEff3P0N_Endcap.SetLineColor(6)
hPtEff3P0N_Endcap.SetLineWidth(2)
hPtEff3P0N_Endcap.SetTitle('3P0N Reconstruction Efficiency vs Pt (Endcap)')
hPtEff3P0N_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N_Endcap.SetStats(0)
hists.append(hPtEff3P0N_Endcap)

hPhiEff3P0N_Endcap = hMatchedTauPhi3P0N_Endcap.Clone('3p0n_tau_phi_eff_endcap')
hPhiEff3P0N_Endcap.Divide(hPhiEff3P0N_Endcap, hTrueTauVisPhi3P0N_Endcap, 1, 1, 'B')
hPhiEff3P0N_Endcap.SetLineColor(418)
hPhiEff3P0N_Endcap.SetLineWidth(2)
hPhiEff3P0N_Endcap.SetTitle('3P0N Reconstruction Efficiency vs Phi (Endcap)')
hPhiEff3P0N_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N_Endcap.SetStats(0)
hists.append(hPhiEff3P0N_Endcap)

hPtEff3P0N_Transition = hMatchedTauPt3P0N_Transition.Clone('3p0n_tau_pt_eff_transition')
hPtEff3P0N_Transition.Divide(hPtEff3P0N_Transition, hTrueTauVisPt3P0N_Transition, 1, 1, 'B')
hPtEff3P0N_Transition.SetLineColor(6)
hPtEff3P0N_Transition.SetLineWidth(2)
hPtEff3P0N_Transition.SetTitle('3P0N Reconstruction Efficiency vs Pt (Transition)')
hPtEff3P0N_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N_Transition.SetStats(0)
hists.append(hPtEff3P0N_Transition)

hPhiEff3P0N_Transition = hMatchedTauPhi3P0N_Transition.Clone('3p0n_tau_phi_eff_transition')
hPhiEff3P0N_Transition.Divide(hPhiEff3P0N_Transition, hTrueTauVisPhi3P0N_Transition, 1, 1, 'B')
hPhiEff3P0N_Transition.SetLineColor(418)
hPhiEff3P0N_Transition.SetLineWidth(2)
hPhiEff3P0N_Transition.SetTitle('3P0N Reconstruction Efficiency vs Phi (Transition)')
hPhiEff3P0N_Transition.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N_Transition.SetStats(0)
hists.append(hPhiEff3P0N_Transition)

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
