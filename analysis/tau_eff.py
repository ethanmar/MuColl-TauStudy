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

def getNRecoPis(reco_tau, relationNavigatorRecoMC):
    n_unique_pis = 0
    n_duplicate_pis = 0
    matched_true_pis = []

    # Loop through associated pfos
    pfos = reco_tau.getParticles()
    for pfo in pfos:

        # Tag charged pions
        if (abs(pfo.getType()) == 211):

            # Get most weighted linked mc particle
            linked_mc_particles = relationNavigatorRecoMC.getRelatedToObjects(pfo)
            linked_weights = relationNavigatorRecoMC.getRelatedToWeights(pfo)
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

                # Sum unique and duplicate matched pions
                if parent_tau is not None:
                    if max_weight_mc_particle not in matched_true_pis:
                        matched_true_pis.append(max_weight_mc_particle)
                        n_unique_pis += 1
                    else:
                        n_duplicate_pis += 1
                    
    return n_unique_pis, n_duplicate_pis

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='tau_eff.root')

args = parser.parse_args()

# Initialize histograms
hists = []

# 1P0N Linked
hLinkedTauPt1P0N = TH1F('linked_1p0n_tau_pt', 'Linked 1P0N Tau Pt', 10, 0, 320)
hists.append(hLinkedTauPt1P0N)

hLinkedTauTheta1P0N = TH1F('linked_1p0n_tau_theta', 'Linked 1P0N Tau Theta', 10, 0, math.pi)
hists.append(hLinkedTauTheta1P0N)

hLinkedTauPhi1P0N = TH1F('linked_1p0n_tau_phi', 'Linked 1P0N Tau Phi', 10, 0, math.pi)
hists.append(hLinkedTauPhi1P0N)

hLinkedTauPt1P0N_Barrel = TH1F('linked_1p0n_tau_pt_barrel', 'Linked 1P0N Tau Pt (Barrel)', 10, 0, 320)
hists.append(hLinkedTauPt1P0N_Barrel)

hLinkedTauPhi1P0N_Barrel = TH1F('linked_1p0n_tau_phi_barrel', 'Linked 1P0N Tau Phi (Barrel)', 10, 0, math.pi)
hists.append(hLinkedTauPhi1P0N_Barrel)

hLinkedTauPt1P0N_CentBarrel = TH1F('linked_1p0n_tau_pt_centbarrel', 'Linked 1P0N Tau Pt (Central Barrel)', 10, 0, 320)
hists.append(hLinkedTauPt1P0N_CentBarrel)

hLinkedTauPhi1P0N_CentBarrel = TH1F('linked_1p0n_tau_phi_centbarrel', 'Linked 1P0N Tau Phi (Central Barrel)', 10, 0, math.pi)
hists.append(hLinkedTauPhi1P0N_CentBarrel)

hLinkedTauPt1P0N_Endcap = TH1F('linked_1p0n_tau_pt_endcap', 'Linked 1P0N Tau Pt (Endcap)', 10, 0, 320)
hists.append(hLinkedTauPt1P0N_Endcap)

hLinkedTauPhi1P0N_Endcap = TH1F('linked_1p0n_tau_phi_endcap', 'Linked 1P0N Tau Phi (Endcap)', 10, 0, math.pi)
hists.append(hLinkedTauPhi1P0N_Endcap)

hLinkedTauPt1P0N_Transition = TH1F('linked_1p0n_tau_pt_transition', 'Linked 1P0N Tau Pt (Transition)', 10, 0, 320)
hists.append(hLinkedTauPt1P0N_Transition)

hLinkedTauPhi1P0N_Transition = TH1F('linked_1p0n_tau_phi_transition', 'Linked 1P0N Tau Phi (Transition)', 10, 0, math.pi)
hists.append(hLinkedTauPhi1P0N_Transition)

# 3P0N Linked
hLinkedTauPt3P0N = TH1F('linked_3p0n_tau_pt', 'Linked 3P0N Tau Pt', 10, 0, 320)
hists.append(hLinkedTauPt3P0N)

hLinkedTauTheta3P0N = TH1F('linked_3p0n_tau_theta', 'Linked 3P0N Tau Theta', 10, 0, math.pi)
hists.append(hLinkedTauTheta3P0N)

hLinkedTauPhi3P0N = TH1F('linked_3p0n_tau_phi', 'Linked 3P0N Tau Phi', 10, 0, math.pi)
hists.append(hLinkedTauPhi3P0N)

hLinkedTauPt3P0N_Barrel = TH1F('linked_3p0n_tau_pt_barrel', 'Linked 3P0N Tau Pt (Barrel)', 10, 0, 320)
hists.append(hLinkedTauPt3P0N_Barrel)

hLinkedTauPhi3P0N_Barrel = TH1F('linked_3p0n_tau_phi_barrel', 'Linked 3P0N Tau Phi (Barrel)', 10, 0, math.pi)
hists.append(hLinkedTauPhi3P0N_Barrel)

hLinkedTauPt3P0N_CentBarrel = TH1F('linked_3p0n_tau_pt_centbarrel', 'Linked 3P0N Tau Pt (Central Barrel)', 10, 0, 320)
hists.append(hLinkedTauPt3P0N_CentBarrel)

hLinkedTauPhi3P0N_CentBarrel = TH1F('linked_3p0n_tau_phi_centbarrel', 'Linked 3P0N Tau Phi (Central Barrel)', 10, 0, math.pi)
hists.append(hLinkedTauPhi3P0N_CentBarrel)

hLinkedTauPt3P0N_Endcap = TH1F('linked_3p0n_tau_pt_endcap', 'Linked 3P0N Tau Pt (Endcap)', 10, 0, 320)
hists.append(hLinkedTauPt3P0N_Endcap)

hLinkedTauPhi3P0N_Endcap = TH1F('linked_3p0n_tau_phi_endcap', 'Linked 3P0N Tau Phi (Endcap)', 10, 0, math.pi)
hists.append(hLinkedTauPhi3P0N_Endcap)

hLinkedTauPt3P0N_Transition = TH1F('linked_3p0n_tau_pt_transition', 'Linked 3P0N Tau Pt (Transition)', 10, 0, 320)
hists.append(hLinkedTauPt3P0N_Transition)

hLinkedTauPhi3P0N_Transition = TH1F('linked_3p0n_tau_phi_transition', 'Linked 3P0N Tau Phi (Transition)', 10, 0, math.pi)
hists.append(hLinkedTauPhi3P0N_Transition)

# 1P0N True
hTrueTauVisPt1P0N = TH1F('true_1p0n_tau_vis_pt', 'True 1P0N Tau Visible Pt', 10, 0, 320)
hists.append(hTrueTauVisPt1P0N)

hTrueTauVisTheta1P0N = TH1F('true_1p0n_tau_vis_theta', 'True 1P0N Tau Visible Theta', 10, 0, math.pi)
hists.append(hTrueTauVisTheta1P0N)

hTrueTauVisPhi1P0N = TH1F('true_1p0n_tau_vis_phi', 'True 1P0N Tau Visible Phi', 10, 0, math.pi)
hists.append(hTrueTauVisPhi1P0N)

hTrueTauVisPt1P0N_Barrel = TH1F('true_1p0n_tau_vis_pt_barrel', 'True 1P0N Tau Visible Pt (Barrel)', 10, 0, 320)
hists.append(hTrueTauVisPt1P0N_Barrel)

hTrueTauVisPhi1P0N_Barrel = TH1F('true_1p0n_tau_vis_phi_barrel', 'True 1P0N Tau Visible Phi (Barrel)', 10, 0, math.pi)
hists.append(hTrueTauVisPhi1P0N_Barrel)

hTrueTauVisPt1P0N_CentBarrel = TH1F('true_1p0n_tau_vis_pt_centbarrel', 'True 1P0N Tau Visible Pt (Central Barrel)', 10, 0, 320)
hists.append(hTrueTauVisPt1P0N_CentBarrel)

hTrueTauVisPhi1P0N_CentBarrel = TH1F('true_1p0n_tau_vis_phi_centbarrel', 'True 1P0N Tau Visible Phi (Central Barrel)', 10, 0, math.pi)
hists.append(hTrueTauVisPhi1P0N_CentBarrel)

hTrueTauVisPt1P0N_Endcap = TH1F('true_1p0n_tau_vis_pt_endcap', 'True 1P0N Tau Visible Pt (Endcap)', 10, 0, 320)
hists.append(hTrueTauVisPt1P0N_Endcap)

hTrueTauVisPhi1P0N_Endcap = TH1F('true_1p0n_tau_vis_phi_endcap', 'True 1P0N Tau Visible Phi (Endcap)', 10, 0, math.pi)
hists.append(hTrueTauVisPhi1P0N_Endcap)

hTrueTauVisPt1P0N_Transition = TH1F('true_1p0n_tau_vis_pt_transition', 'True 1P0N Tau Visible Pt (Transition)', 10, 0, 320)
hists.append(hTrueTauVisPt1P0N_Transition)

hTrueTauVisPhi1P0N_Transition = TH1F('true_1p0n_tau_vis_phi_transition', 'True 1P0N Tau Visible Phi (Transition)', 10, 0, math.pi)
hists.append(hTrueTauVisPhi1P0N_Transition)

# 3P0N True
hTrueTauVisPt3P0N = TH1F('true_3p0n_tau_vis_pt', 'True 3P0N Tau Visible Pt', 10, 0, 320)
hists.append(hTrueTauVisPt3P0N)

hTrueTauVisTheta3P0N = TH1F('true_3p0n_tau_vis_theta', 'True 3P0N Tau Visible Theta', 10, 0, math.pi)
hists.append(hTrueTauVisTheta3P0N)

hTrueTauVisPhi3P0N = TH1F('true_3p0n_tau_vis_phi', 'True 3P0N Tau Visible Phi', 10, 0, math.pi)
hists.append(hTrueTauVisPhi3P0N)

hTrueTauVisPt3P0N_Barrel = TH1F('true_3p0n_tau_vis_pt_barrel', 'True 3P0N Tau Visible Pt (Barrel)', 10, 0, 320)
hists.append(hTrueTauVisPt3P0N_Barrel)

hTrueTauVisPhi3P0N_Barrel = TH1F('true_3p0n_tau_vis_phi_barrel', 'True 3P0N Tau Visible Phi (Barrel)', 10, 0, math.pi)
hists.append(hTrueTauVisPhi3P0N_Barrel)

hTrueTauVisPt3P0N_CentBarrel = TH1F('true_3p0n_tau_vis_pt_centbarrel', 'True 3P0N Tau Visible Pt (Central Barrel)', 10, 0, 320)
hists.append(hTrueTauVisPt3P0N_CentBarrel)

hTrueTauVisPhi3P0N_CentBarrel = TH1F('true_3p0n_tau_vis_phi_centbarrel', 'True 3P0N Tau Visible Phi (Central Barrel)', 10, 0, math.pi)
hists.append(hTrueTauVisPhi3P0N_CentBarrel)

hTrueTauVisPt3P0N_Endcap = TH1F('true_3p0n_tau_vis_pt_endcap', 'True 3P0N Tau Visible Pt (Endcap)', 10, 0, 320)
hists.append(hTrueTauVisPt3P0N_Endcap)

hTrueTauVisPhi3P0N_Endcap = TH1F('true_3p0n_tau_vis_phi_endcap', 'True 3P0N Tau Visible Phi (Endcap)', 10, 0, math.pi)
hists.append(hTrueTauVisPhi3P0N_Endcap)

hTrueTauVisPt3P0N_Transition = TH1F('true_3p0n_tau_vis_pt_transition', 'True 3P0N Tau Visible Pt (Transition)', 10, 0, 320)
hists.append(hTrueTauVisPt3P0N_Transition)

hTrueTauVisPhi3P0N_Transition = TH1F('true_3p0n_tau_vis_phi_transition', 'True 3P0N Tau Visible Phi (Transition)', 10, 0, math.pi)
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

n_pi_true = []
n_pi_reco = []

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
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauPFOLink')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        # Instantiate relation navigators to parse tauReco and RecoMC links
        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        # Loop through mc particles
        for mcParticle in mcParticles:

            # Tag taus
            if (abs(mcParticle.getPDG()) == 15):
                
                # Get decay mode
                decay_mode = getDecayMode(mcParticle)

                # Get visible properties
                vis_props = getVisibleProperties(mcParticle)
                vis_pt = getPt(vis_props)
                vis_theta = getTheta(vis_props)
                vis_phi = getPhi(vis_props)
                
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
        
        # Loop through reco taus
        for tau in taus:

            # Get linked true tau
            true_tau = getLinkedMCTau(tau, relationNavigatorTau, relationNavigatorRecoMC)

            # Get true decay mode
            decay_mode = getDecayMode(true_tau)

            # Get visible properties of linked true tau
            vis_props = getVisibleProperties(true_tau)
            vis_pt = getPt(vis_props)
            vis_theta = getTheta(vis_props)
            vis_phi = getPhi(vis_props)

            # Get number of associated pions
            n_unique_pis, n_duplicate_pis = getNRecoPis(tau, relationNavigatorRecoMC)

            if (decay_mode == 0 or decay_mode == 4):
                n_pi_reco.append(n_unique_pis + n_duplicate_pis)
                if (decay_mode == 0):
                    n_pi_true.append(1)
                elif (decay_mode == 4):
                    n_pi_true.append(3)
            
            # Fill linked 1P0N hists
            if (n_unique_pis == 1 and n_duplicate_pis == 0 and decay_mode == 0):
                n_1p0n_reco += 1
                hLinkedTauPt1P0N.Fill(vis_pt)
                hLinkedTauTheta1P0N.Fill(vis_theta)
                hLinkedTauPhi1P0N.Fill(vis_phi)

                if (vis_theta > 0.70 and vis_theta < 2.45):
                    hLinkedTauPt1P0N_Barrel.Fill(vis_pt)
                    hLinkedTauPhi1P0N_Barrel.Fill(vis_phi)

                if (vis_theta > 1 and vis_theta < 2):
                    hLinkedTauPt1P0N_CentBarrel.Fill(vis_pt)
                    hLinkedTauPhi1P0N_CentBarrel.Fill(vis_phi)

                elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                    hLinkedTauPt1P0N_Transition.Fill(vis_pt)
                    hLinkedTauPhi1P0N_Transition.Fill(vis_phi)

                elif (vis_theta < 0.577 or vis_theta > 2.56):
                    hLinkedTauPt1P0N_Endcap.Fill(vis_pt)
                    hLinkedTauPhi1P0N_Endcap.Fill(vis_phi)

            # Fill linked 3P0N hists
            elif (n_unique_pis == 3 and n_duplicate_pis == 0 and decay_mode == 4):
                n_3p0n_reco += 1
                hLinkedTauPt3P0N.Fill(vis_pt)
                hLinkedTauTheta3P0N.Fill(vis_theta)
                hLinkedTauPhi3P0N.Fill(vis_phi)
                
                if (vis_theta > 0.70 and vis_theta < 2.45):
                    hLinkedTauPt3P0N_Barrel.Fill(vis_pt)
                    hLinkedTauPhi3P0N_Barrel.Fill(vis_phi)

                if (vis_theta > 1 and vis_theta < 2):
                    hLinkedTauPt3P0N_CentBarrel.Fill(vis_pt)
                    hLinkedTauPhi3P0N_CentBarrel.Fill(vis_phi)

                elif ((vis_theta > 0.577 and vis_theta < 1.0) or (vis_theta > 2.0 and vis_theta < 2.56)):
                    hLinkedTauPt3P0N_Transition.Fill(vis_pt)
                    hLinkedTauPhi3P0N_Transition.Fill(vis_phi)

                elif (vis_theta < 0.577 or vis_theta > 2.56):
                    hLinkedTauPt3P0N_Endcap.Fill(vis_pt)
                    hLinkedTauPhi3P0N_Endcap.Fill(vis_phi)

    reader.close()

print(f'1p0n eff: {n_1p0n_reco/n_1p0n_true}')
print(f'3p0n eff: {n_3p0n_reco/n_3p0n_true}')
    
# Create confusion matrices
plt.figure()

plt.clf()
disp = ConfusionMatrixDisplay.from_predictions(n_pi_reco, n_pi_true, normalize='pred', include_values=True, cmap='Blues', colorbar=True)
disp.ax_.set_ylabel(r'Number of Reconstructed Prongs (0 $\pi^0$s)', fontsize=12)
disp.ax_.set_xlabel(r'Number of True Prongs (0 $\pi^0$s)', fontsize=12)
disp.ax_.invert_yaxis()
plt.tight_layout()
plt.savefig('tau_confusion.png')

# Create 1P0N efficiency hists
hPtEff1P0N = hLinkedTauPt1P0N.Clone('1p0n_tau_pt_eff')
hPtEff1P0N.Divide(hPtEff1P0N, hTrueTauVisPt1P0N, 1, 1, 'B')
hPtEff1P0N.SetLineColor(6)
hPtEff1P0N.SetLineWidth(2)
hPtEff1P0N.SetTitle('1P0N Reconstruction Efficiency vs Pt')
hPtEff1P0N.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N.SetStats(0)
hists.append(hPtEff1P0N)

hThetaEff1P0N = hLinkedTauTheta1P0N.Clone('1p0n_tau_theta_eff')
hThetaEff1P0N.Divide(hThetaEff1P0N, hTrueTauVisTheta1P0N, 1, 1, 'B')
hThetaEff1P0N.SetLineColor(7)
hThetaEff1P0N.SetLineWidth(2)
hThetaEff1P0N.SetTitle('1P0N Reconstruction Efficiency vs Theta')
hThetaEff1P0N.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hThetaEff1P0N.GetYaxis().SetTitle('#epsilon')
hThetaEff1P0N.SetStats(0)
hists.append(hThetaEff1P0N)

hPhiEff1P0N = hLinkedTauPhi1P0N.Clone('1p0n_tau_phi_eff')
hPhiEff1P0N.Divide(hPhiEff1P0N, hTrueTauVisPhi1P0N, 1, 1, 'B')
hPhiEff1P0N.SetLineColor(418)
hPhiEff1P0N.SetLineWidth(2)
hPhiEff1P0N.SetTitle('1P0N Reconstruction Efficiency vs Phi')
hPhiEff1P0N.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N.SetStats(0)
hists.append(hPhiEff1P0N)

hPtEff1P0N_Barrel = hLinkedTauPt1P0N_Barrel.Clone('1p0n_tau_pt_eff_barrel')
hPtEff1P0N_Barrel.Divide(hPtEff1P0N_Barrel, hTrueTauVisPt1P0N_Barrel, 1, 1, 'B')
hPtEff1P0N_Barrel.SetLineColor(6)
hPtEff1P0N_Barrel.SetLineWidth(2)
hPtEff1P0N_Barrel.SetTitle('1P0N Reconstruction Efficiency vs Pt (Barrel)')
hPtEff1P0N_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N_Barrel.SetStats(0)
hists.append(hPtEff1P0N_Barrel)

hPhiEff1P0N_Barrel = hLinkedTauPhi1P0N_Barrel.Clone('1p0n_tau_phi_eff_barrel')
hPhiEff1P0N_Barrel.Divide(hPhiEff1P0N_Barrel, hTrueTauVisPhi1P0N_Barrel, 1, 1, 'B')
hPhiEff1P0N_Barrel.SetLineColor(418)
hPhiEff1P0N_Barrel.SetLineWidth(2)
hPhiEff1P0N_Barrel.SetTitle('1P0N Reconstruction Efficiency vs Phi (Barrel)')
hPhiEff1P0N_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N_Barrel.SetStats(0)
hists.append(hPhiEff1P0N_Barrel)

hPtEff1P0N_CentBarrel = hLinkedTauPt1P0N_CentBarrel.Clone('1p0n_tau_pt_eff_centbarrel')
hPtEff1P0N_CentBarrel.Divide(hPtEff1P0N_CentBarrel, hTrueTauVisPt1P0N_CentBarrel, 1, 1, 'B')
hPtEff1P0N_CentBarrel.SetLineColor(6)
hPtEff1P0N_CentBarrel.SetLineWidth(2)
hPtEff1P0N_CentBarrel.SetTitle('1P0N Reconstruction Efficiency vs Pt (Central Barrel)')
hPtEff1P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N_CentBarrel.SetStats(0)
hists.append(hPtEff1P0N_CentBarrel)

hPhiEff1P0N_CentBarrel = hLinkedTauPhi1P0N_CentBarrel.Clone('1p0n_tau_phi_eff_centbarrel')
hPhiEff1P0N_CentBarrel.Divide(hPhiEff1P0N_CentBarrel, hTrueTauVisPhi1P0N_CentBarrel, 1, 1, 'B')
hPhiEff1P0N_CentBarrel.SetLineColor(418)
hPhiEff1P0N_CentBarrel.SetLineWidth(2)
hPhiEff1P0N_CentBarrel.SetTitle('1P0N Reconstruction Efficiency vs Phi (Central Barrel)')
hPhiEff1P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N_CentBarrel.SetStats(0)
hists.append(hPhiEff1P0N_CentBarrel)

hPtEff1P0N_Endcap = hLinkedTauPt1P0N_Endcap.Clone('1p0n_tau_pt_eff_endcap')
hPtEff1P0N_Endcap.Divide(hPtEff1P0N_Endcap, hTrueTauVisPt1P0N_Endcap, 1, 1, 'B')
hPtEff1P0N_Endcap.SetLineColor(6)
hPtEff1P0N_Endcap.SetLineWidth(2)
hPtEff1P0N_Endcap.SetTitle('1P0N Reconstruction Efficiency vs Pt (Endcap)')
hPtEff1P0N_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N_Endcap.SetStats(0)
hists.append(hPtEff1P0N_Endcap)

hPhiEff1P0N_Endcap = hLinkedTauPhi1P0N_Endcap.Clone('1p0n_tau_phi_eff_endcap')
hPhiEff1P0N_Endcap.Divide(hPhiEff1P0N_Endcap, hTrueTauVisPhi1P0N_Endcap, 1, 1, 'B')
hPhiEff1P0N_Endcap.SetLineColor(418)
hPhiEff1P0N_Endcap.SetLineWidth(2)
hPhiEff1P0N_Endcap.SetTitle('1P0N Reconstruction Efficiency vs Phi (Endcap)')
hPhiEff1P0N_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N_Endcap.SetStats(0)
hists.append(hPhiEff1P0N_Endcap)

hPtEff1P0N_Transition = hLinkedTauPt1P0N_Transition.Clone('1p0n_tau_pt_eff_transition')
hPtEff1P0N_Transition.Divide(hPtEff1P0N_Transition, hTrueTauVisPt1P0N_Transition, 1, 1, 'B')
hPtEff1P0N_Transition.SetLineColor(6)
hPtEff1P0N_Transition.SetLineWidth(2)
hPtEff1P0N_Transition.SetTitle('1P0N Reconstruction Efficiency vs Pt (Transition)')
hPtEff1P0N_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff1P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPtEff1P0N_Transition.SetStats(0)
hists.append(hPtEff1P0N_Transition)

hPhiEff1P0N_Transition = hLinkedTauPhi1P0N_Transition.Clone('1p0n_tau_phi_eff_transition')
hPhiEff1P0N_Transition.Divide(hPhiEff1P0N_Transition, hTrueTauVisPhi1P0N_Transition, 1, 1, 'B')
hPhiEff1P0N_Transition.SetLineColor(418)
hPhiEff1P0N_Transition.SetLineWidth(2)
hPhiEff1P0N_Transition.SetTitle('1P0N Reconstruction Efficiency vs Phi (Transition)')
hPhiEff1P0N_Transition.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff1P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPhiEff1P0N_Transition.SetStats(0)
hists.append(hPhiEff1P0N_Transition)

# Create reco-three-prong efficiency hists
hPtEff3P0N = hLinkedTauPt3P0N.Clone('3p0n_tau_pt_eff')
hPtEff3P0N.Divide(hPtEff3P0N, hTrueTauVisPt3P0N, 1, 1, 'B')
hPtEff3P0N.SetLineColor(6)
hPtEff3P0N.SetLineWidth(2)
hPtEff3P0N.SetTitle('3P0N Reconstruction Efficiency vs Pt')
hPtEff3P0N.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N.SetStats(0)
hists.append(hPtEff3P0N)

hThetaEff3P0N = hLinkedTauTheta3P0N.Clone('3p0n_tau_theta_eff')
hThetaEff3P0N.Divide(hThetaEff3P0N, hTrueTauVisTheta3P0N, 1, 1, 'B')
hThetaEff3P0N.SetLineColor(7)
hThetaEff3P0N.SetLineWidth(2)
hThetaEff3P0N.SetTitle('3P0N Reconstruction Efficiency vs Theta')
hThetaEff3P0N.GetXaxis().SetTitle('True Visible Tau #theta [rad]')
hThetaEff3P0N.GetYaxis().SetTitle('#epsilon')
hThetaEff3P0N.SetStats(0)
hists.append(hThetaEff3P0N)

hPhiEff3P0N = hLinkedTauPhi3P0N.Clone('3p0n_tau_phi_eff')
hPhiEff3P0N.Divide(hPhiEff3P0N, hTrueTauVisPhi3P0N, 1, 1, 'B')
hPhiEff3P0N.SetLineColor(418)
hPhiEff3P0N.SetLineWidth(2)
hPhiEff3P0N.SetTitle('3P0N Reconstruction Efficiency vs Phi')
hPhiEff3P0N.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N.SetStats(0)
hists.append(hPhiEff3P0N)

hPtEff3P0N_Barrel = hLinkedTauPt3P0N_Barrel.Clone('3p0n_tau_pt_eff_barrel')
hPtEff3P0N_Barrel.Divide(hPtEff3P0N_Barrel, hTrueTauVisPt3P0N_Barrel, 1, 1, 'B')
hPtEff3P0N_Barrel.SetLineColor(6)
hPtEff3P0N_Barrel.SetLineWidth(2)
hPtEff3P0N_Barrel.SetTitle('3P0N Reconstruction Efficiency vs Pt (Barrel)')
hPtEff3P0N_Barrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N_Barrel.SetStats(0)
hists.append(hPtEff3P0N_Barrel)

hPhiEff3P0N_Barrel = hLinkedTauPhi3P0N_Barrel.Clone('3p0n_tau_phi_eff_barrel')
hPhiEff3P0N_Barrel.Divide(hPhiEff3P0N_Barrel, hTrueTauVisPhi3P0N_Barrel, 1, 1, 'B')
hPhiEff3P0N_Barrel.SetLineColor(418)
hPhiEff3P0N_Barrel.SetLineWidth(2)
hPhiEff3P0N_Barrel.SetTitle('3P0N Reconstruction Efficiency vs Phi (Barrel)')
hPhiEff3P0N_Barrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N_Barrel.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N_Barrel.SetStats(0)
hists.append(hPhiEff3P0N_Barrel)

hPtEff3P0N_CentBarrel = hLinkedTauPt3P0N_CentBarrel.Clone('3p0n_tau_pt_eff_centbarrel')
hPtEff3P0N_CentBarrel.Divide(hPtEff3P0N_CentBarrel, hTrueTauVisPt3P0N_CentBarrel, 1, 1, 'B')
hPtEff3P0N_CentBarrel.SetLineColor(6)
hPtEff3P0N_CentBarrel.SetLineWidth(2)
hPtEff3P0N_CentBarrel.SetTitle('3P0N Reconstruction Efficiency vs Pt (Central Barrel)')
hPtEff3P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N_CentBarrel.SetStats(0)
hists.append(hPtEff3P0N_CentBarrel)

hPhiEff3P0N_CentBarrel = hLinkedTauPhi3P0N_CentBarrel.Clone('3p0n_tau_phi_eff_centbarrel')
hPhiEff3P0N_CentBarrel.Divide(hPhiEff3P0N_CentBarrel, hTrueTauVisPhi3P0N_CentBarrel, 1, 1, 'B')
hPhiEff3P0N_CentBarrel.SetLineColor(418)
hPhiEff3P0N_CentBarrel.SetLineWidth(2)
hPhiEff3P0N_CentBarrel.SetTitle('3P0N Reconstruction Efficiency vs Phi (Central Barrel)')
hPhiEff3P0N_CentBarrel.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N_CentBarrel.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N_CentBarrel.SetStats(0)
hists.append(hPhiEff3P0N_CentBarrel)

hPtEff3P0N_Endcap = hLinkedTauPt3P0N_Endcap.Clone('3p0n_tau_pt_eff_endcap')
hPtEff3P0N_Endcap.Divide(hPtEff3P0N_Endcap, hTrueTauVisPt3P0N_Endcap, 1, 1, 'B')
hPtEff3P0N_Endcap.SetLineColor(6)
hPtEff3P0N_Endcap.SetLineWidth(2)
hPtEff3P0N_Endcap.SetTitle('3P0N Reconstruction Efficiency vs Pt (Endcap)')
hPtEff3P0N_Endcap.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N_Endcap.SetStats(0)
hists.append(hPtEff3P0N_Endcap)

hPhiEff3P0N_Endcap = hLinkedTauPhi3P0N_Endcap.Clone('3p0n_tau_phi_eff_endcap')
hPhiEff3P0N_Endcap.Divide(hPhiEff3P0N_Endcap, hTrueTauVisPhi3P0N_Endcap, 1, 1, 'B')
hPhiEff3P0N_Endcap.SetLineColor(418)
hPhiEff3P0N_Endcap.SetLineWidth(2)
hPhiEff3P0N_Endcap.SetTitle('3P0N Reconstruction Efficiency vs Phi (Endcap)')
hPhiEff3P0N_Endcap.GetXaxis().SetTitle('True Visible Tau #phi [rad]')
hPhiEff3P0N_Endcap.GetYaxis().SetTitle('#epsilon')
hPhiEff3P0N_Endcap.SetStats(0)
hists.append(hPhiEff3P0N_Endcap)

hPtEff3P0N_Transition = hLinkedTauPt3P0N_Transition.Clone('3p0n_tau_pt_eff_transition')
hPtEff3P0N_Transition.Divide(hPtEff3P0N_Transition, hTrueTauVisPt3P0N_Transition, 1, 1, 'B')
hPtEff3P0N_Transition.SetLineColor(6)
hPtEff3P0N_Transition.SetLineWidth(2)
hPtEff3P0N_Transition.SetTitle('3P0N Reconstruction Efficiency vs Pt (Transition)')
hPtEff3P0N_Transition.GetXaxis().SetTitle('True Visible Tau Pt [GeV/c]')
hPtEff3P0N_Transition.GetYaxis().SetTitle('#epsilon')
hPtEff3P0N_Transition.SetStats(0)
hists.append(hPtEff3P0N_Transition)

hPhiEff3P0N_Transition = hLinkedTauPhi3P0N_Transition.Clone('3p0n_tau_phi_eff_transition')
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
