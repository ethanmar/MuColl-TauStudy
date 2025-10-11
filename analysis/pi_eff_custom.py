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

from tau_mc_link import getDecayMode

def getVisibleProperties(mc_tau):

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

def isTauDaughter(mc_pi):
    parents = mc_pi.getParents()
    for parent in parents:
        if (abs(parent.getPDG()) == 15):
            return True
    return False

def getRecoPis(pfos, mc_particles):
    # unique_matched_mc_pis = []
    # duplicate_matched_mc_pis = []
    matched_mc_pis = []
    matched_reco_pis = []
    
    # Loop over pfos
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
                
                # Compute energy resolution
                mc_energy = deltaR_mc_pi.getEnergy()
                energy_dif = abs(pfo_energy - mc_energy)

                # Mc pi matched if it has smallest energy resolution
                if (energy_dif < min_energy_dif):
                    min_energy_dif = energy_dif
                    matched_mc_pi = deltaR_mc_pi

            if (matched_mc_pi is not None):
                matched_mc_pis.append(matched_mc_pi)
                matched_reco_pis.append(pfo)
                '''if (matched_mc_pi not in unique_matched_mc_pis):
                    unique_matched_mc_pis.append(matched_mc_pi)
                else:
                    duplicate_matched_mc_pis.append(matched_mc_pi)'''

    return matched_mc_pis, matched_reco_pis

def getRecoOthers(pfos, mc_particles):
    matched_mc_others = []
    matched_reco_others = []
    
    # Loop over pfos
    for pfo in pfos:

        # Tag not charged pions
        if (abs(pfo.getType()) != 211):

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
            deltaR_mc_others = []
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
                        deltaR_mc_others.append(mc_particle)
                        
            # Loop over mc particles with deltaR less than 0.1 rad
            min_energy_dif = 1e6
            matched_mc_other = None
            for deltaR_mc_other in deltaR_mc_others:
                
                # Compute energy resolution
                mc_energy = deltaR_mc_other.getEnergy()
                energy_dif = abs(pfo_energy - mc_energy)

                # Mc pi matched if it has smallest energy resolution
                if (energy_dif < min_energy_dif):
                    min_energy_dif = energy_dif
                    matched_mc_other = deltaR_mc_other

            if (matched_mc_other is not None):
                matched_mc_others.append(matched_mc_other)
                matched_reco_others.append(pfo)

    return matched_mc_others, matched_reco_others

def getPtBin(pt, pt_bins):
    for i in range(len(pt_bins)-1):
        if (pt >= pt_bins[i] and pt < pt_bins[i+1]):
            return i
    return len(pt_bins)-2

def getThetaBin(theta, theta_bins):
    for i in range(len(theta_bins)-1):
        if (theta >= theta_bins[i] and theta < theta_bins[i+1]):
            return i
    return len(theta_bins)-2

def getMinPtBin(min_pt, min_pt_bins):
    for i in range(len(min_pt_bins)-1):
        if (min_pt >= min_pt_bins[i] and min_pt < min_pt_bins[i+1]):
            return i
    return len(min_pt_bins)-2

def getAvgAngle(mc_tau):
    pis = []
    daughters = mc_tau.getDaughters()
    for daughter in daughters:
        if (abs(daughter.getPDG()) == 211):
            pis.append(daughter)

    angle_sum = 0
    total = 0
    for i in range(0, len(pis)-1):
        imom = pis[i].getMomentum()
        ip = np.array([imom[0], imom[1], imom[2]])
        inorm = np.linalg.norm(ip)
        for j in range(1+i, len(pis)):
            jmom = pis[j].getMomentum()
            jp = np.array([jmom[0], jmom[1], jmom[2]])
            jnorm = np.linalg.norm(jp)
            dot = np.dot(ip, jp)
            cos_angle = dot/(inorm*jnorm)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            angle = np.arccos(cos_angle)
            angle_sum += angle
            total += 1
            
    return angle_sum/total

def getMinPt(mc_tau):
    min_pt = 1e6
    daughters = mc_tau.getDaughters()
    for daughter in daughters:
        if (abs(daughter.getPDG()) == 211):
            mom = daughter.getMomentum()
            pt = math.sqrt(mom[0]**2 + mom[1]**2)
            if (pt < min_pt):
                min_pt = pt
    return min_pt

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='pi_eff_custom.root')

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

hPiAvgAngleTrue3P0N = TH1F('3p0n_pi_true_avg_angle', 'True 3P0N Pion Average Angle', 20, 0, 0.30)
hists.append(hPiAvgAngleTrue3P0N)

hPiPtMatched3P0N = TH1F('3p0n_pi_matched_pT', 'Matched 3P0N Pion Pt', 10, 0, 320)
hists.append(hPiPtMatched3P0N)

hPiPhiMatched3P0N = TH1F('3p0n_pi_matched_phi', 'Matched 3P0N Pion Phi', 10, 0, math.pi)
hists.append(hPiPhiMatched3P0N)

hPiThetaMatched3P0N = TH1F('3p0n_pi_matched_theta', 'Matched 3P0N Pion Theta', 10, 0, math.pi)
hists.append(hPiThetaMatched3P0N)

hPiAvgAngleMatched3P0N = TH1F('3p0n_pi_matched_avg_angle', 'Matched 3P0N Pion Average Angle', 20, 0, 0.30)
hists.append(hPiAvgAngleMatched3P0N)

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

# Failed 3P0N Pions
hFailedPiPtTrue3P0N = TH1F('3p0n_failed_pi_true_pt', 'True 3P0N Failed Pion Pt', 50, 0, 320)
hists.append(hFailedPiPtTrue3P0N)

hFailedPiThetaTrue3P0N = TH1F('3p0n_failed_pi_true_theta', 'True 3P0N Failed Pion Theta', 50, 0, math.pi/2)
hists.append(hFailedPiThetaTrue3P0N)

hFailedPiPhiTrue3P0N = TH1F('3p0n_failed_pi_true_phi', 'True 3P0N Failed Pion Phi', 50, 0, math.pi)
hists.append(hFailedPiPhiTrue3P0N)

hFailedPiMaxPtRatioTrue3P0N = TH1F('3p0n_failed_pi_true_max_pt_ratio', 'True 3P0N Failed Pion Max Pt Ratio', 50, 0, 1)
hists.append(hFailedPiMaxPtRatioTrue3P0N)

hFailedPiAvgAngleTrue3P0N = TH1F('3p0n_failed_pi_true_avg_angle', 'True 3P0N Failed Pion Average Angle', 50, 0, 0.25)
hists.append(hFailedPiAvgAngleTrue3P0N)

# Passed 3P0N Pions
hPassedPiPtTrue3P0N = TH1F('3p0n_passed_pi_true_pt', 'True 3P0N Passed Pion Pt', 50, 0, 320)
hists.append(hPassedPiPtTrue3P0N)

hPassedPiThetaTrue3P0N = TH1F('3p0n_passed_pi_true_theta', 'True 3P0N Passed Pion Theta', 50, 0, math.pi/2)
hists.append(hPassedPiThetaTrue3P0N)

hPassedPiPhiTrue3P0N = TH1F('3p0n_passed_pi_true_phi', 'True 3P0N Passed Pion Phi', 50, 0, math.pi)
hists.append(hPassedPiPhiTrue3P0N)

hPassedPiMaxPtRatioTrue3P0N = TH1F('3p0n_passed_pi_true_max_pt_ratio', 'True 3P0N Passed Pion Max Pt Ratio', 50, 0, 1)
hists.append(hPassedPiMaxPtRatioTrue3P0N)

hPassedPiAvgAngleTrue3P0N = TH1F('3p0n_passed_pi_true_avg_angle', 'True 3P0N Passed Pion Average Angle', 50, 0, 0.25)
hists.append(hPassedPiAvgAngleTrue3P0N)

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

# n_pi_true = []
# n_pi_reco = []
pi_cm = np.zeros((5, 3))
# n_pi_unique = []
# n_pi_duplicate = []

theta_bins = np.linspace(0, np.pi/2, 11)
pt_bins = np.linspace(0, 320, 11)
min_pt_bins = np.linspace(0, 50, 11)
numerator_1p0n = np.zeros((len(pt_bins)-1, len(theta_bins)-1))
denominator_1p0n = np.zeros((len(pt_bins)-1, len(theta_bins)-1))
numerator_3p0n = np.zeros((len(pt_bins)-1, len(theta_bins)-1))
denominator_3p0n = np.zeros((len(pt_bins)-1, len(theta_bins)-1))
numerator_3p0n_min_pt = np.zeros((len(min_pt_bins)-1, len(theta_bins)-1))
denominator_3p0n_min_pt = np.zeros((len(min_pt_bins)-1, len(theta_bins)-1))

n_electrons_1p0n = np.zeros((len(pt_bins)-1, len(theta_bins)-1))
n_electrons_3p0n = np.zeros((len(pt_bins)-1, len(theta_bins)-1))
counts_1p0n = np.zeros((len(pt_bins)-1, len(theta_bins)-1))
counts_3p0n = np.zeros((len(pt_bins)-1, len(theta_bins)-1))

n_electrons_1p0n_min_pt = np.zeros((len(min_pt_bins)-1, len(theta_bins)-1))
n_electrons_3p0n_min_pt = np.zeros((len(min_pt_bins)-1, len(theta_bins)-1))
counts_1p0n_min_pt = np.zeros((len(min_pt_bins)-1, len(theta_bins)-1))
counts_3p0n_min_pt = np.zeros((len(min_pt_bins)-1, len(theta_bins)-1))

event_file = "failed_events.txt"

# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        pfos = event.getCollection('PandoraPFOs')
        mc_particles = event.getCollection('MCParticle')

        # True decay mode
        decay_mode = -1
        
        # Visible observables
        vis_props = []
        vis_pt = 0
        vis_theta = 0
        vis_phi = 0

        # Average angle
        avg_angle = 0

        # True min pt of 3p0n
        min_pt = 1e6
        
        # Loop through mc_particles
        for mc_particle in mc_particles:

            # Tag taus
            if (abs(mc_particle.getPDG()) == 15):

                # Get decay mode
                decay_mode = getDecayMode(mc_particle)
                
                # Get visible properties
                vis_props = getVisibleProperties(mc_particle)
                vis_pt = getPt(vis_props)
                vis_theta = getTheta(vis_props)
                vis_phi = getPhi(vis_props)

                if decay_mode == 4:
                    avg_angle = getAvgAngle(mc_particle)

                min_pt = getMinPt(mc_particle)

        pt_bin = getPtBin(vis_pt, pt_bins)
        theta_bin = getThetaBin(abs(vis_theta-(math.pi/2)), theta_bins)
        min_pt_bin = getMinPtBin(min_pt, min_pt_bins)
        
        if (decay_mode == 0):
            denominator_1p0n[pt_bin][theta_bin] += 1
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

        elif (decay_mode == 4):
            denominator_3p0n[pt_bin][theta_bin] += 1
            if min_pt <= 50:
                denominator_3p0n_min_pt[min_pt_bin][theta_bin] += 1
            n_3p0n_true += 1
            hPiPtTrue3P0N.Fill(vis_pt)
            hPiThetaTrue3P0N.Fill(vis_theta)
            hPiPhiTrue3P0N.Fill(vis_phi)
            hPiAvgAngleTrue3P0N.Fill(avg_angle)

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

        # Get number of matched pions
        matched_mc_pis, matched_reco_pis = getRecoPis(pfos, mc_particles)
        n_matched_pis = len(matched_mc_pis)

        # Get other matched pfos
        matched_mc_others, matched_reco_others = getRecoOthers(pfos, mc_particles)

        n_electrons = 0
        for matched_reco_other in matched_reco_others:
            if (abs(matched_reco_other.getType()) == 11):
                n_electrons += 1

        if decay_mode == 0:
            n_electrons_1p0n[pt_bin][theta_bin] += n_electrons
            counts_1p0n[pt_bin][theta_bin] += 1
            if min_pt <= 50:
                n_electrons_1p0n_min_pt[min_pt_bin][theta_bin] += n_electrons
                counts_1p0n_min_pt[min_pt_bin][theta_bin] += 1
        elif decay_mode == 4:
            n_electrons_3p0n[pt_bin][theta_bin] += n_electrons
            counts_3p0n[pt_bin][theta_bin] += 1
            if min_pt <= 50:
                n_electrons_3p0n_min_pt[min_pt_bin][theta_bin] += n_electrons
                counts_3p0n_min_pt[min_pt_bin][theta_bin] += 1
        
        if (decay_mode == 0 or decay_mode == 4):
            if n_matched_pis <= 4:
                if decay_mode == 0:
                    pi_cm[n_matched_pis][0] += 1
                elif decay_mode == 4:
                    pi_cm[n_matched_pis][2] += 1
            # n_pi_unique.append(n_unique_pis)
            # n_pi_duplicate.append(n_duplicate_pis)
                
        if (n_matched_pis == 1 and decay_mode == 0):
            numerator_1p0n[pt_bin][theta_bin] += 1
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

        elif (n_matched_pis == 3 and decay_mode == 4):
            numerator_3p0n[pt_bin][theta_bin] += 1
            if min_pt <= 50:
                numerator_3p0n_min_pt[min_pt_bin][theta_bin] += 1
            n_3p0n_reco += 1
            hPiPtMatched3P0N.Fill(vis_pt)
            hPiThetaMatched3P0N.Fill(vis_theta)
            hPiPhiMatched3P0N.Fill(vis_phi)
            hPiAvgAngleMatched3P0N.Fill(avg_angle)

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

            pt_ratio_max = 0
            for mc_particle in mc_particles:

                if (abs(mc_particle.getPDG()) == 211 and isTauDaughter(mc_particle)):
                    mom = mc_particle.getMomentum()
                    px = mom[0]
                    py = mom[1]
                    pz = mom[2]
                    pt = math.sqrt(px**2 + py**2)
                    theta =  math.atan2(pt, pz)
                    phi = math.atan2(py, px)
                    hPassedPiPtTrue3P0N.Fill(pt)
                    hPassedPiThetaTrue3P0N.Fill(abs(theta-(math.pi/2)))
                    hPassedPiPhiTrue3P0N.Fill(phi)

                    pt_ratio = pt/vis_pt
                    if pt_ratio > pt_ratio_max:
                        pt_ratio_max = pt_ratio
            hPassedPiMaxPtRatioTrue3P0N.Fill(pt_ratio_max)
            hPassedPiAvgAngleTrue3P0N.Fill(avg_angle)

            '''if (vis_pt > 0 and vis_pt < 64):
                if (vis_theta > 1.0 and vis_theta < 2.0):
                    with open(event_file, "a") as f:
                        f.write(f"{event.getEventNumber()}\n")'''

        elif (n_matched_pis != 3 and decay_mode == 4):

            pt_ratio_max = 0
            for mc_particle in mc_particles:

                if (abs(mc_particle.getPDG()) == 211 and isTauDaughter(mc_particle)):
                    mom = mc_particle.getMomentum()
                    px = mom[0]
                    py = mom[1]
                    pz = mom[2]
                    pt = math.sqrt(px**2 + py**2)
                    theta =  math.atan2(pt, pz)
                    phi = math.atan2(py, px)
                    hFailedPiPtTrue3P0N.Fill(pt)
                    hFailedPiThetaTrue3P0N.Fill(abs(theta-(math.pi/2)))
                    hFailedPiPhiTrue3P0N.Fill(phi)

                    pt_ratio = pt/vis_pt
                    if pt_ratio > pt_ratio_max:
                        pt_ratio_max = pt_ratio
            hFailedPiMaxPtRatioTrue3P0N.Fill(pt_ratio_max)
            hFailedPiAvgAngleTrue3P0N.Fill(avg_angle)
            
            if (vis_pt > 0 and vis_pt < 64):
                if ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                    with open(event_file, "a") as f:
                        f.write(f"{event.getEventNumber()}\n")
                        
    reader.close()

print(f'1p0n eff: {n_1p0n_reco/n_1p0n_true}')
print(f'3p0n eff: {n_3p0n_reco/n_3p0n_true}')

# Create confusion matrices
plt.figure()

pi_column_sums = pi_cm.sum(axis=0, keepdims=True)
pi_cm_norm = np.divide(pi_cm.astype('float'), pi_column_sums, out=np.zeros_like(pi_cm, dtype=float), where=pi_column_sums != 0)

plt.imshow(pi_cm_norm, origin='lower', cmap='Blues')
plt.colorbar()
plt.xlabel(r'Number of True $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
plt.ylabel(r'Number of Reconstructed $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
plt.xticks(np.arange(pi_cm_norm.shape[1]), labels=['1', '2', '3'])
plt.yticks(np.arange(pi_cm_norm.shape[0]), labels=['0', '1', '2', '3', '4'])
for i in range(pi_cm_norm.shape[0]):
    for j in range(pi_cm_norm.shape[1]):
        value = pi_cm_norm[i,j]
        plt.text(j, i, f"{value:.2f}", ha="center", va="center", color="black")
plt.tight_layout()
plt.savefig('pi_confusion.png')

print(f'1p0n_denominator: {denominator_1p0n}')
print(f'\n3p0n_denominator: {denominator_3p0n}')

# efficiency_1p0n = numerator_1p0n/denominator_1p0n
# efficiency_3p0n = numerator_3p0n/denominator_3p0n
efficiency_1p0n = np.divide(numerator_1p0n, denominator_1p0n, out=np.zeros_like(numerator_1p0n, dtype=float), where=denominator_1p0n != 0)
efficiency_3p0n = np.divide(numerator_3p0n, denominator_3p0n, out=np.zeros_like(numerator_3p0n, dtype=float), where=denominator_3p0n != 0)
# error_1p0n = np.sqrt((efficiency_1p0n*(1-efficiency_1p0n))/denominator_1p0n)
# error_3p0n = np.sqrt((efficiency_3p0n*(1-efficiency_3p0n))/denominator_3p0n)
numer_1p0n = efficiency_1p0n*(1-efficiency_1p0n)
error_1p0n = np.sqrt(np.divide(numer_1p0n, denominator_1p0n, out=np.zeros_like(numer_1p0n, dtype=float), where=denominator_1p0n != 0))
numer_3p0n = efficiency_3p0n*(1-efficiency_3p0n)
error_3p0n = np.sqrt(np.divide(numer_3p0n, denominator_3p0n, out=np.zeros_like(numer_3p0n, dtype=float), where=denominator_3p0n != 0))

pt_centers = 0.5*(pt_bins[:-1] + pt_bins[1:])
theta_centers = 0.5*(theta_bins[:-1] + theta_bins[1:])

plt.clf()
plt.imshow(efficiency_1p0n, extent=[theta_bins[0], theta_bins[-1], pt_bins[0], pt_bins[-1]], aspect='auto', origin='lower', cmap='viridis')
plt.xticks(theta_bins, labels=[f"{t:.2f}" for t in theta_bins])
plt.yticks(pt_bins, labels=[f"{pt:.2f}" for pt in pt_bins])
plt.xlabel("Visible Tau Theta [rad]")
plt.ylabel("Visible Tau Pt [GeV/c]")
plt.colorbar(label=" 1P0N Pion Efficiency")
n_rows, n_columns = efficiency_1p0n.shape
for i in range(n_rows):
    for j in range(n_columns):
        x = theta_centers[j]
        y = pt_centers[i]
        value = efficiency_1p0n[i][j]
        err = error_1p0n[i][j]
        plt.text(x, y, rf"${value:.2f}$" + "\n" + rf"$\pm {err:.2f}$", ha='center', va='center', fontsize=6, color='black')
plt.tight_layout()
plt.savefig('1p0n_eff_3d_hist_high_bins.png')

plt.clf()
plt.imshow(efficiency_3p0n, extent=[theta_bins[0], theta_bins[-1], pt_bins[0], pt_bins[-1]], aspect='auto', origin='lower', cmap='viridis')
plt.xticks(theta_bins, labels=[f"{t:.2f}" for t in theta_bins])
plt.yticks(pt_bins, labels=[f"{pt:.2f}" for pt in pt_bins])
plt.xlabel("Visible Tau Theta [rad]")
plt.ylabel("Visible Tau Pt [GeV/c]")
plt.colorbar(label=" 3P0N Pion Efficiency")
n_rows, n_columns = efficiency_3p0n.shape
for i in range(n_rows):
    for j in range(n_columns):
        x = theta_centers[j]
        y = pt_centers[i]
        value = efficiency_3p0n[i][j]
        err = error_3p0n[i][j]
        plt.text(x, y, rf"${value:.2f}$" + "\n" + rf"$\pm {err:.2f}$", ha='center', va='center', fontsize=6, color='black')
plt.tight_layout()
plt.savefig('3p0n_eff_3d_hist_high_bins.png')

efficiency_3p0n_min_pt = np.divide(numerator_3p0n_min_pt, denominator_3p0n_min_pt, out=np.zeros_like(numerator_3p0n_min_pt, dtype=float), where=denominator_3p0n_min_pt != 0)
numer_3p0n_min_pt = efficiency_3p0n_min_pt*(1-efficiency_3p0n_min_pt)
error_3p0n_min_pt = np.sqrt(np.divide(numer_3p0n_min_pt, denominator_3p0n_min_pt, out=np.zeros_like(numer_3p0n_min_pt, dtype=float), where=denominator_3p0n_min_pt != 0))
min_pt_centers = 0.5*(min_pt_bins[:-1] + min_pt_bins[1:])

plt.clf()
plt.imshow(efficiency_3p0n_min_pt, extent=[theta_bins[0], theta_bins[-1], min_pt_bins[0], min_pt_bins[-1]], aspect='auto', origin='lower', cmap='viridis')
plt.xticks(theta_bins, labels=[f"{t:.2f}" for t in theta_bins])
plt.yticks(min_pt_bins, labels=[f"{m:.2f}" for m in min_pt_bins])
plt.xlabel("Visible Tau Theta [rad]")
plt.ylabel("Visible Tau Min Pt [GeV/c]")
plt.colorbar(label=" 3P0N Pion Efficiency")
n_rows, n_columns = efficiency_3p0n_min_pt.shape
for i in range(n_rows):
    for j in range(n_columns):
        x = theta_centers[j]
        y = min_pt_centers[i]
        value = efficiency_3p0n_min_pt[i][j]
        err = error_3p0n_min_pt[i][j]
        plt.text(x, y, rf"${value:.2f}$" + "\n" + rf"$\pm {err:.2f}$", ha='center', va='center', fontsize=6, color='black')
plt.tight_layout()
plt.savefig('3p0n_min_pt_eff_3d_hist_high_bins.png')

n_electrons_1p0n_norm = np.divide(n_electrons_1p0n, counts_1p0n, out=np.zeros_like(n_electrons_1p0n, dtype=float), where=counts_1p0n != 0)
n_electrons_3p0n_norm = np.divide(n_electrons_3p0n, counts_3p0n, out=np.zeros_like(n_electrons_3p0n, dtype=float), where=counts_3p0n != 0)

plt.clf()
plt.imshow(n_electrons_1p0n_norm, extent=[theta_bins[0], theta_bins[-1], pt_bins[0], pt_bins[-1]], aspect='auto', origin='lower', cmap='viridis')
plt.xticks(theta_bins, labels=[f"{t:.2f}" for t in theta_bins])
plt.yticks(pt_bins, labels=[f"{pt:.2f}" for pt in pt_bins])
plt.xlabel(r"Visible True 1P0N Tau $|\theta - \pi/2|$ [rad]")
plt.ylabel(r"Visible True 1P0N Tau $p_T$ [GeV/c]")
plt.colorbar(label="Normalized Number of Reco Electrons")
n_rows, n_columns = n_electrons_1p0n_norm.shape
for i in range(n_rows):
    for j in range(n_columns):
        x = theta_centers[j]
        y = pt_centers[i]
        value = n_electrons_1p0n_norm[i][j]
        plt.text(x, y, rf"${value:.2f}$", ha='center', va='center', fontsize=6, color='black')
plt.tight_layout()
plt.savefig('1p0n_n_electrons_norm.png')

plt.clf()
plt.imshow(n_electrons_3p0n_norm, extent=[theta_bins[0], theta_bins[-1], pt_bins[0], pt_bins[-1]], aspect='auto', origin='lower', cmap='viridis')
plt.xticks(theta_bins, labels=[f"{t:.2f}" for t in theta_bins])
plt.yticks(pt_bins, labels=[f"{pt:.2f}" for pt in pt_bins])
plt.xlabel(r"Visible True 3P0N Tau $|\theta - \pi/2|$ [rad]")
plt.ylabel(r"Visible True 3P0N Tau $p_T$ [GeV/c]")
plt.colorbar(label="Normalized Number of Reco Electrons")
n_rows, n_columns = n_electrons_3p0n_norm.shape
for i in range(n_rows):
    for j in range(n_columns):
        x = theta_centers[j]
        y = pt_centers[i]
        value = n_electrons_3p0n_norm[i][j]
        plt.text(x, y, rf"${value:.2f}$", ha='center', va='center', fontsize=6, color='black')
plt.tight_layout()
plt.savefig('3p0n_n_electrons_norm.png')

n_electrons_1p0n_min_pt_norm = np.divide(n_electrons_1p0n_min_pt, counts_1p0n_min_pt, out=np.zeros_like(n_electrons_1p0n_min_pt, dtype=float), where=counts_1p0n_min_pt != 0)
n_electrons_3p0n_min_pt_norm = np.divide(n_electrons_3p0n_min_pt, counts_3p0n_min_pt, out=np.zeros_like(n_electrons_3p0n_min_pt, dtype=float), where=counts_3p0n_min_pt != 0)

plt.clf()
plt.imshow(n_electrons_1p0n_min_pt_norm, extent=[theta_bins[0], theta_bins[-1], min_pt_bins[0], min_pt_bins[-1]], aspect='auto', origin='lower', cmap='viridis')
plt.xticks(theta_bins, labels=[f"{t:.2f}" for t in theta_bins])
plt.yticks(min_pt_bins, labels=[f"{pt:.2f}" for pt in min_pt_bins])
plt.xlabel("Visible Tau Theta [rad]")
plt.ylabel("Visible Tau Min Pt [GeV/c]")
plt.colorbar(label="Normalized Number of Electrons (1P0N)")
n_rows, n_columns = n_electrons_1p0n_min_pt_norm.shape
for i in range(n_rows):
    for j in range(n_columns):
        x = theta_centers[j]
        y = min_pt_centers[i]
        value = n_electrons_1p0n_min_pt_norm[i][j]
        plt.text(x, y, rf"${value:.2f}$", ha='center', va='center', fontsize=6, color='black')
plt.tight_layout()
plt.savefig('1p0n_n_electrons_min_pt_norm.png')

plt.clf()
plt.imshow(n_electrons_3p0n_min_pt_norm, extent=[theta_bins[0], theta_bins[-1], min_pt_bins[0], min_pt_bins[-1]], aspect='auto', origin='lower', cmap='viridis')
plt.xticks(theta_bins, labels=[f"{t:.2f}" for t in theta_bins])
plt.yticks(min_pt_bins, labels=[f"{pt:.2f}" for pt in min_pt_bins])
plt.xlabel("Visible Tau Theta [rad]")
plt.ylabel("Visible Tau Min Pt [GeV/c]")
plt.colorbar(label="Normalized Number of Electrons (3P0N)")
n_rows, n_columns = n_electrons_3p0n_min_pt_norm.shape
for i in range(n_rows):
    for j in range(n_columns):
        x = theta_centers[j]
        y = min_pt_centers[i]
        value = n_electrons_3p0n_min_pt_norm[i][j]
        plt.text(x, y, rf"${value:.2f}$", ha='center', va='center', fontsize=6, color='black')
plt.tight_layout()
plt.savefig('3p0n_n_electrons_min_pt_norm.png')

'''plt.clf()
disp = ConfusionMatrixDisplay.from_predictions(n_pi_unique, n_pi_true, normalize='pred', include_values=True, cmap='Blues', colorbar=True)
disp.ax_.set_ylabel(r'Number of Unique Reco $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
disp.ax_.set_xlabel(r'Number of True $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
disp.ax_.invert_yaxis()
plt.tight_layout()
plt.savefig('pi_confusion_unique.png')

plt.clf()
disp = ConfusionMatrixDisplay.from_predictions(n_pi_duplicate, n_pi_true, normalize='pred', include_values=True, cmap='Blues', colorbar=True)
disp.ax_.set_ylabel(r'Number of Duplicate Reco $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
disp.ax_.set_xlabel(r'Number of True $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
disp.ax_.invert_yaxis()
plt.tight_layout()
plt.savefig('pi_confusion_duplicate.png')'''

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

hPiAvgAngleEff3P0N = hPiAvgAngleMatched3P0N.Clone('3p0n_pi_avg_angle_eff')
hPiAvgAngleEff3P0N.Divide(hPiAvgAngleEff3P0N, hPiAvgAngleTrue3P0N, 1, 1, 'B')
hPiAvgAngleEff3P0N.SetLineColor(438)
hPiAvgAngleEff3P0N.SetLineWidth(2)
hPiAvgAngleEff3P0N.SetTitle('')
hPiAvgAngleEff3P0N.GetXaxis().SetTitle('True Average Angle [rad]')
hPiAvgAngleEff3P0N.GetYaxis().SetTitle('#epsilon')
hPiAvgAngleEff3P0N.SetStats(0)
hists.append(hPiAvgAngleEff3P0N)

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
