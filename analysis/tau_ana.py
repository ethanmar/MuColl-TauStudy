from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np

from tau_mc_link import getLinkedMCTau, getVisibleProperties, getDecayMode, getNRecoQPis

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='tau_ana.root')

args = parser.parse_args()

# Initialize histograms
hists = []

# Tau hists

# Reco
hTauE = TH1F('tau_energy', 'Reconstructed Tau Energy', 500, 0, 1500)
hists.append(hTauE)

hTauPt = TH1F('tau_pt', 'Reconstructed Tau Pt', 20, 0, 350)
hists.append(hTauPt)

hTauTheta = TH1F('tau_theta', 'Reconstructed Tau Theta', 20, 0, math.pi)
hists.append(hTauTheta)

hTauPhi = TH1F('tau_phi', 'Reconstructed Tau Phi', 20, 0, math.pi)
hists.append(hTauPhi)

hTauNDaughters = TH1F('n_tau_daughters', 'Number of Reconstructed Tau Daughters', 10, 0, 10)
hists.append(hTauNDaughters)

hTauDaughterType = TH1F('tau_daughter_types', 'Reconstructed Tau Daughter Types', 1100, 0, 2200)
hists.append(hTauDaughterType)

hLinkedMCTauPt = TH1F('linked_mc_tau_pt', 'Linked MC Tau Pt', 20, 0, 320)
hists.append(hLinkedMCTauPt)

hLinkedMCTauTheta = TH1F('linked_mc_tau_theta', 'Linked MC Tau Theta', 20, 0, math.pi)
hists.append(hLinkedMCTauTheta)

hLinkedMCTauPhi = TH1F('linked_mc_tau_phi', 'Linked MC Tau Phi', 20, 0, math.pi)
hists.append(hLinkedMCTauPhi)

# One-prong
hLinkedMCTauPt1P = TH1F('linked_1p_mc_tau_pt', 'Linked One-Prong MC Tau Pt', 20, 0, 320)
hists.append(hLinkedMCTauPt1P)

hLinkedMCTauTheta1P = TH1F('linked_1p_mc_tau_theta', 'Linked One-Prong MC Tau Theta', 20, 0, math.pi)
hists.append(hLinkedMCTauTheta1P)

hLinkedMCTauPhi1P = TH1F('linked_1p_mc_tau_phi', 'Linked One-Prong MC Tau Phi', 20, 0, math.pi)
hists.append(hLinkedMCTauPhi1P)

# Three-prong
hLinkedMCTauPt3P = TH1F('linked_3p_mc_tau_pt', 'Linked Three-Prong MC Tau Pt', 20, 0, 320)
hists.append(hLinkedMCTauPt3P)

hLinkedMCTauTheta3P = TH1F('linked_3p_mc_tau_theta', 'Linked Three-Prong MC Tau Theta', 20, 0, math.pi)
hists.append(hLinkedMCTauTheta3P)

hLinkedMCTauPhi3P = TH1F('linked_3p_mc_tau_phi', 'Linked Three-Prong MC Tau Phi', 20, 0, math.pi)
hists.append(hLinkedMCTauPhi3P)

# One-prong (reco)
hLinkedMCTauPt1PReco = TH1F('linked_1p_reco_mc_tau_pt', 'Linked Reco-One-Prong MC Tau Pt', 20, 0, 320)
hists.append(hLinkedMCTauPt1PReco)

hLinkedMCTauTheta1PReco = TH1F('linked_1p_reco_mc_tau_theta', 'Linked Reco-One-Prong MC Tau Theta', 20, 0, math.pi)
hists.append(hLinkedMCTauTheta1PReco)

hLinkedMCTauPhi1PReco = TH1F('linked_1p_reco_mc_tau_phi', 'Linked Reco-One-Prong MC Tau Phi', 20, 0, math.pi)
hists.append(hLinkedMCTauPhi1PReco)

hTauE1P = TH1F('1p_tau_energy', 'Reconstructed One-Prong Tau Energy', 500, 0, 1500)
hists.append(hTauE1P)

hTauPt1P = TH1F('1p_tau_pt', 'Reconstructed One-Prong Tau Pt', 20, 0, 350)
hists.append(hTauPt1P)

hTauTheta1P = TH1F('1p_tau_theta', 'Reconstructed One-Prong Tau Theta', 20, 0, math.pi)
hists.append(hTauTheta1P)

hTauPhi1P = TH1F('1p_tau_phi', 'Reconstructed One-Prong Tau Phi', 20, 0, math.pi)
hists.append(hTauPhi1P)

hTauNDaughters1P = TH1F('n_1p_tau_daughters', 'Number of Reconstructed One-Prong Tau Daughters', 10, 0, 10)
hists.append(hTauNDaughters1P)

hTauDaughterType1P = TH1F('1p_tau_daughter_types', 'Reconstructed One-Prong Tau Daughter Types', 2200, 0, 2200)
hists.append(hTauDaughterType1P)

hTauResE1P = TH1F('1p_tau_res_E', 'One-Prong Tau Residual Energy', 60, -60, 60)
hTauResE1P.GetXaxis().SetTitle(r'E_{reco} - E_{true} [GeV]')
hists.append(hTauResE1P)

hTauResPt1P = TH1F('1p_tau_res_pt', 'One-Prong Tau Residual Pt', 60, -60, 60)
hTauResPt1P.GetXaxis().SetTitle(r'Pt_{reco} - Pt_{true} [GeV/c]')
hists.append(hTauResPt1P)

# Three-prong (reco)
hLinkedMCTauPt3PReco = TH1F('linked_3p_reco_mc_tau_pt', 'Linked Reco-Three-Prong MC Tau Pt', 20, 0, 320)
hists.append(hLinkedMCTauPt3PReco)

hLinkedMCTauTheta3PReco = TH1F('linked_3p_reco_mc_tau_theta', 'Linked Reco-Three-Prong MC Tau Theta', 20, 0, math.pi)
hists.append(hLinkedMCTauTheta3PReco)

hLinkedMCTauPhi3PReco = TH1F('linked_3p_reco_mc_tau_phi', 'Linked Reco-Three-Prong MC Tau Phi', 20, 0, math.pi)
hists.append(hLinkedMCTauPhi3PReco)

hTauE3P = TH1F('3p_tau_energy', 'Reconstructed Three-Prong Tau Energy', 500, 0, 1500)
hists.append(hTauE3P)

hTauPt3P = TH1F('3p_tau_pt', 'Reconstructed Three-Prong Tau Pt', 20, 0, 350)
hists.append(hTauPt3P)

hTauTheta3P = TH1F('3p_tau_theta', 'Reconstructed Three-Prong Tau Theta', 20, 0, math.pi)
hists.append(hTauTheta3P)

hTauPhi3P = TH1F('3p_tau_phi', 'Reconstructed Three-Prong Tau Phi', 20, 0, math.pi)
hists.append(hTauPhi3P)

hTauNDaughters3P = TH1F('3p_n_tau_daughters', 'Number of Reconstructed Three-Prong Tau Daughters', 10, 0, 10)
hists.append(hTauNDaughters3P)

hTauDaughterType3P = TH1F('3p_tau_daughter_types', 'Reconstructed Three-Prong Tau Daughter Types', 2200, 0, 2200)
hists.append(hTauDaughterType3P)

hTauResE3P = TH1F('3p_tau_res_E', 'Three-Prong Tau Residual Energy', 200, -100, 100)
hists.append(hTauResE3P)

hTauResPt3P = TH1F('3p_tau_res_pT', 'Three-Prong Tau Residual Pt', 200, -100, 100)
hists.append(hTauResPt3P)

#Truth
hTauVisETrue = TH1F('tau_true_energy', 'True Tau Visible Energy', 500, 0, 1500)
hists.append(hTauVisETrue)

hTauVisPtTrue = TH1F('tau_true_pT', 'True Tau Visible Pt', 20, 0, 320)
hists.append(hTauVisPtTrue)

hTauVisThetaTrue = TH1F('tau_true_theta', 'True Tau Visible Theta', 20, 0, math.pi)
hists.append(hTauVisThetaTrue)

hTauVisPhiTrue = TH1F('tau_phi_true', 'True Tau Visible Phi', 20, 0, math.pi)
hists.append(hTauVisPhiTrue)

hTauNVisDaughtersTrue = TH1F('n_tau_visible_daughters_true', 'Number of True Visible Tau Daughters', 10, 0, 10)
hists.append(hTauNVisDaughtersTrue)

hTauVisDaughterTypeTrue = TH1F('visible_tau_daughter_types_true', 'True Visible Tau Daughter Types', 300, 0, 300)
hists.append(hTauVisDaughterTypeTrue)

# One-prong
hTauVisPtTrue1P = TH1F('1p_true_vis_pT', 'True One-Prong Tau Visible Pt', 20, 0, 320)
hists.append(hTauVisPtTrue1P)

hTauVisThetaTrue1P = TH1F('1p_true_vis_theta', 'True One-Prong Tau Visible Theta', 20, 0, math.pi)
hists.append(hTauVisThetaTrue1P)

hTauVisPhiTrue1P = TH1F('1p_true_vis_phi', 'True One-Prong Tau Visible Phi', 20, 0, math.pi)
hists.append(hTauVisPhiTrue1P)

hTauTotPtTrue1P = TH1F('1p_true_tot_pT', 'True One-Prong Tau Total Pt', 20, 0, 320)
hists.append(hTauTotPtTrue1P)

hTauTotThetaTrue1P = TH1F('1p_true_tot_theta', 'True One-Prong Tau Total Theta', 20, 0, math.pi)
hists.append(hTauTotThetaTrue1P)

hTauTotPhiTrue1P = TH1F('1p_true_tot_phi', 'True One-Prong Tau Total Phi', 20, 0, math.pi)
hists.append(hTauTotPhiTrue1P)

# Three-prong
hTauVisPtTrue3P = TH1F('3p_tau_true_pT', 'True Three-Prong Tau Visible Pt', 20, 0, 320)
hists.append(hTauVisPtTrue3P)

hTauVisThetaTrue3P = TH1F('3p_tau_true_theta', 'True Three-Prong Tau Visible Theta', 20, 0, math.pi)
hists.append(hTauVisThetaTrue3P)

hTauVisPhiTrue3P = TH1F('3p_tau_phi_true', 'True Three-Prong Tau Visible Phi', 20, 0, math.pi)
hists.append(hTauVisPhiTrue3P)

hTauTotPtTrue3P = TH1F('3p_tau_tot_pT', 'True Three-Prong Tau Total Pt', 20, 0, 320)
hists.append(hTauTotPtTrue3P)

hTauTotThetaTrue3P = TH1F('3p_tau_tot_theta', 'True Three-Prong Tau Total Theta', 20, 0, math.pi)
hists.append(hTauTotThetaTrue3P)

hTauTotPhiTrue3P = TH1F('3p_tau_phi_tot', 'True Three-Prong Tau Total Phi', 20, 0, math.pi)
hists.append(hTauTotPhiTrue3P)

# 1P Pions
hPiPtTrue1P = TH1F('1p_pi_true_pT', 'True One-Prong Pion Pt', 20, 0, 320)
hists.append(hPiPtTrue1P)

hPiPhiTrue1P = TH1F('1p_pi_true_phi', 'True One-Prong Pion Phi', 20, 0, math.pi)
hists.append(hPiPhiTrue1P)

hPiThetaTrue1P = TH1F('1p_pi_true_theta', 'True One-Prong Pion Theta', 20, 0, math.pi)
hists.append(hPiThetaTrue1P)

hPiPtMatched1P = TH1F('1p_pi_matched_pT', 'Matched One-Prong Pion Pt', 20, 0, 320)
hists.append(hPiPtMatched1P)

hPiPhiMatched1P = TH1F('1p_pi_matched_phi', 'Matched One-Prong Pion Phi', 20, 0, math.pi)
hists.append(hPiPhiMatched1P)

hPiThetaMatched1P = TH1F('1p_pi_matched_theta', 'Matched One-Prong Pion Theta', 20, 0, math.pi)
hists.append(hPiThetaMatched1P)

# 3P Pions
hPiPtTrue3P = TH1F('3p_pi_true_pT', 'True Three-Prong Pion Pt', 20, 0, 320)
hists.append(hPiPtTrue3P)

hPiPhiTrue3P = TH1F('3p_pi_true_phi', 'True Three-Prong Pion Phi', 20, 0, math.pi)
hists.append(hPiPhiTrue3P)

hPiThetaTrue3P = TH1F('3p_pi_true_theta', 'True Three-Prong Pion Theta', 20, 0, math.pi)
hists.append(hPiThetaTrue3P)

hPiPtMatched3P = TH1F('3p_pi_matched_pT', 'Matched Three-Prong Pion Pt', 20, 0, 320)
hists.append(hPiPtMatched3P)

hPiPhiMatched3P = TH1F('3p_pi_matched_phi', 'Matched Three-Prong Pion Phi', 20, 0, math.pi)
hists.append(hPiPhiMatched3P)

hPiThetaMatched3P = TH1F('3p_pi_matched_theta', 'Matched Three-Prong Pion Theta', 20, 0, math.pi)
hists.append(hPiThetaMatched3P)

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

# Keep track of one-prong reco daughter types
one_prong_types = {}
    
# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        taus = event.getCollection('TauRec_PFO')
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauRecLink_PFO')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        # Instantiate relation navigators to parse tauReco and RecoMC links
        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)
        
        # Loop through tau PFOs
        for tau in taus:
            
            # Instantiate and calculate tau observables
            E = tau.getEnergy()
            p = tau.getMomentum()
            px = p[0]
            py = p[1]
            pz = p[2]
            pt = math.sqrt(px**2 + py**2)
            theta = math.acos(p[2]/(math.sqrt(pt**2+p[2]**2)))
            phi = math.acos(p[0]/pt)
            daughters = tau.getParticles()

            # Fill tau hists
            hTauE.Fill(E)
            hTauPt.Fill(pt)
            hTauTheta.Fill(theta)
            hTauPhi.Fill(phi)
            hTauNDaughters.Fill(len(daughters))

            # Loop over reco daughters and fill reco daughter hists
            for daughter in daughters:
                type_ = abs(daughter.getType())
                hTauDaughterType.Fill(type_)

            
            # Get linked MC tau
            mcTau = getLinkedMCTau(tau, relationNavigatorTau, relationNavigatorRecoMC)

            # Get visible properties of linked MC tau
            E_vis, px_vis, py_vis, pz_vis, n_daughters_vis, vis_daughter_types = getVisibleProperties(mcTau)
            pt_vis = math.sqrt(px_vis**2 + py_vis**2)
            theta_vis = math.acos(pz_vis/(math.sqrt(pt_vis**2 + pz_vis**2)))
            phi_vis = math.acos(px_vis/pt_vis)

            # Fill linked MC tau hists
            hLinkedMCTauPt.Fill(pt_vis)
            hLinkedMCTauTheta.Fill(theta_vis)
            hLinkedMCTauPhi.Fill(phi_vis)

            # Get linked MC tau decay mode
            decayMode = getDecayMode(mcTau)

            # Fill one-prong hists
            if (decayMode == 0):
                hLinkedMCTauPt1P.Fill(pt_vis)
                hLinkedMCTauTheta1P.Fill(theta_vis)
                hLinkedMCTauPhi1P.Fill(phi_vis)

            # Fill three-prong hists
            elif (decayMode == 4):
                hLinkedMCTauPt3P.Fill(pt_vis)
                hLinkedMCTauTheta3P.Fill(theta_vis)
                hLinkedMCTauPhi3P.Fill(phi_vis)

            # Get number of reco charged pions in reco tau daughters
            nRecoQPis = getNRecoQPis(tau)

            # Fill correct reco-one-prong hists
            if (nRecoQPis == 1 and decayMode == 0):
                hLinkedMCTauPt1PReco.Fill(pt_vis)
                hLinkedMCTauTheta1PReco.Fill(theta_vis)
                hLinkedMCTauPhi1PReco.Fill(phi_vis)                

            # Fill reco-one-prong hists
            if (nRecoQPis == 1 and decayMode == 0):
                hTauE1P.Fill(E)
                hTauPt1P.Fill(pt)
                hTauTheta1P.Fill(theta)
                hTauPhi1P.Fill(phi)
                hTauNDaughters1P.Fill(len(daughters))
                for daughter in daughters:
                    type_ = abs(daughter.getType())
                    hTauDaughterType1P.Fill(type_)
                    if str(type_) in one_prong_types:
                        one_prong_types[str(type_)] += 1
                    else:
                        one_prong_types[str(type_)] = 1

            if (nRecoQPis == 1 and decayMode == 0):
                hTauResE1P.Fill(E-E_vis)
                hTauResPt1P.Fill(pt-pt_vis)

            # Fill correct reco-three-prong hists
            if (nRecoQPis == 3 and decayMode == 4):
                hLinkedMCTauPt3PReco.Fill(pt_vis)
                hLinkedMCTauTheta3PReco.Fill(theta_vis)
                hLinkedMCTauPhi3PReco.Fill(phi_vis)

            # Fill reco-three-prong hists
            if (nRecoQPis == 3):
                hTauE3P.Fill(E)
                hTauPt3P.Fill(pt)
                hTauTheta3P.Fill(theta)
                hTauPhi3P.Fill(phi)
                hTauNDaughters3P.Fill(len(daughters))
                for daughter in daughters:
                    hTauDaughterType3P.Fill(abs(daughter.getType()))
                hTauResE3P.Fill(E-E_vis)
                hTauResPt3P.Fill(pt-pt_vis)
                    
        # Loop through MC particles
        for mcParticle in mcParticles:

            # Get MC taus
            pdg = abs(mcParticle.getPDG())
            if pdg == 15:

                # Get visible properties
                E_vis, px_vis, py_vis, pz_vis, n_daughters_vis, vis_daughter_types = getVisibleProperties(mcParticle)
                pt_vis = math.sqrt(px_vis**2 + py_vis**2)
                theta_vis = math.acos(pz_vis/(math.sqrt(pt_vis**2 + pz_vis**2)))
                phi_vis = math.acos(px_vis/pt_vis)

                # Fill visible mc tau hists
                hTauVisETrue.Fill(E_vis)
                hTauVisPtTrue.Fill(pt_vis)
                hTauVisThetaTrue.Fill(theta_vis)
                hTauVisPhiTrue.Fill(phi_vis)
                hTauNVisDaughtersTrue.Fill(n_daughters_vis)
                for vis_daughter_type in vis_daughter_types:
                    hTauVisDaughterTypeTrue.Fill(vis_daughter_type)

                # Get total properties
                tot_mom = mcParticle.getMomentum()
                tot_px = tot_mom[0]
                tot_py = tot_mom[1]
                tot_pz = tot_mom[2]
                tot_pt = math.sqrt(tot_px**2 + tot_py**2)
                tot_theta = math.acos(tot_pz/(math.sqrt(tot_pt**2 + tot_pz**2)))
                tot_phi = math.acos(tot_px/tot_pt)

                # Get mc particle decay mode
                decayMode = getDecayMode(mcParticle)

                # Fill one-prong hists
                if (decayMode == 0):
                    hTauVisPtTrue1P.Fill(pt_vis)
                    hTauVisThetaTrue1P.Fill(theta_vis)
                    hTauVisPhiTrue1P.Fill(phi_vis)
                    hTauTotPtTrue1P.Fill(tot_pt)
                    hTauTotThetaTrue1P.Fill(tot_theta)
                    hTauTotPhiTrue1P.Fill(tot_phi)

                # Fill three-prong hists
                elif (decayMode == 4):
                    hTauVisPtTrue3P.Fill(pt_vis)
                    hTauVisThetaTrue3P.Fill(theta_vis)
                    hTauVisPhiTrue3P.Fill(phi_vis)
                    hTauTotPtTrue3P.Fill(tot_pt)
                    hTauTotThetaTrue3P.Fill(tot_theta)
                    hTauTotPhiTrue3P.Fill(tot_phi)

                # Investigate pion reco for 1P decays
                if (decayMode == 0):

                    # Loop over mc 1P tau daughters
                    daughters = mcParticle.getDaughters()
                    for daughter in daughters:
                        if abs(daughter.getPDG()) == 211:
                            hPiPtTrue1P.Fill(pt_vis)
                            hPiPhiTrue1P.Fill(phi_vis)
                            hPiThetaTrue1P.Fill(theta_vis)

                    # Loop over pandora pfos
                    matched_mc_pi = None
                    for pfo in pfos:
                        # Obtain linked mc pions
                        if abs(pfo.getType()) == 211:
                            linkedMCParticles = relationNavigatorRecoMC.getRelatedToObjects(pfo)
                            for linkedMCParticle in linkedMCParticles:
                                if abs(linkedMCParticle.getPDG()) == 211:
                                    matched_mc_pi = linkedMCParticle

                    if matched_mc_pi is not None:
                        hPiPtMatched1P.Fill(pt_vis)
                        hPiPhiMatched1P.Fill(phi_vis)
                        hPiThetaMatched1P.Fill(theta_vis)

                # Investigate pion reco for 3P decays
                if (decayMode == 4):
                    hPiPtTrue3P.Fill(pt_vis)
                    hPiPhiTrue3P.Fill(phi_vis)
                    hPiThetaTrue3P.Fill(theta_vis)

                    # Loop over pandora pfos
                    matched_mc_pi = 0
                    for pfo in pfos:
                        # Obtain linked mc pions
                        if abs(pfo.getType()) == 211:
                            linkedMCParticles = relationNavigatorRecoMC.getRelatedToObjects(pfo)
                            for linkedMCParticle in linkedMCParticles:
                                if abs(linkedMCParticle.getPDG()) == 211:
                                    matched_mc_pi += 1

                    if matched_mc_pi == 3:
                        hPiPtMatched3P.Fill(pt_vis)
                        hPiPhiMatched3P.Fill(phi_vis)
                        hPiThetaMatched3P.Fill(theta_vis)

    # Close file
    reader.close()

# Create one-prong efficiency hists
hPtEff1P = hLinkedMCTauPt1P.Clone('1p_pt_eff')
hPtEff1P.Divide(hPtEff1P, hTauVisPtTrue1P, 1, 1, 'B')
hPtEff1P.SetLineColor(6)
hPtEff1P.SetLineWidth(2)
hPtEff1P.SetTitle('One-Prong Reconstruction Efficiency vs Pt')
hPtEff1P.GetXaxis().SetTitle('True Visible Pt [GeV/c]')
hPtEff1P.GetYaxis().SetTitle('#epsilon')
hPtEff1P.SetStats(0)
hists.append(hPtEff1P)

hThetaEff1P = hLinkedMCTauTheta1P.Clone('1p_theta_eff')
hThetaEff1P.Divide(hThetaEff1P, hTauVisThetaTrue1P, 1, 1, 'B')
hThetaEff1P.SetLineColor(7)
hThetaEff1P.SetLineWidth(2)
hThetaEff1P.SetTitle('One-Prong Reconstruction Efficiency vs Theta')
hThetaEff1P.GetXaxis().SetTitle('True Visible #theta [rad]')
hThetaEff1P.GetYaxis().SetTitle('#epsilon')
hThetaEff1P.SetStats(0)
hists.append(hThetaEff1P)

hPhiEff1P = hLinkedMCTauPhi1P.Clone('1p_phi_eff')
hPhiEff1P.Divide(hPhiEff1P, hTauVisPhiTrue1P, 1, 1, 'B')
hPhiEff1P.SetLineColor(418)
hPhiEff1P.SetLineWidth(2)
hPhiEff1P.SetTitle('One-Prong Reconstruction Efficiency vs Phi')
hPhiEff1P.GetXaxis().SetTitle('True Visible #phi [rad]')
hPhiEff1P.GetYaxis().SetTitle('#epsilon')
hPhiEff1P.SetStats(0)
hists.append(hPhiEff1P)

# Create three-prong efficiency hists
hPtEff3P = hLinkedMCTauPt3P.Clone('3p_pt_eff')
hPtEff3P.Divide(hPtEff3P, hTauVisPtTrue3P, 1, 1, 'B')
hPtEff3P.SetLineColor(6)
hPtEff3P.SetLineWidth(2)
hPtEff3P.SetTitle('Three-Prong Reconstruction Efficiency vs Pt')
hPtEff3P.GetXaxis().SetTitle('True Visible Pt [GeV/c]')
hPtEff3P.GetYaxis().SetTitle('#epsilon')
hPtEff3P.SetStats(0)
hists.append(hPtEff3P)

hThetaEff3P = hLinkedMCTauTheta3P.Clone('3p_theta_eff')
hThetaEff3P.Divide(hThetaEff3P, hTauVisThetaTrue3P, 1, 1, 'B')
hThetaEff3P.SetLineColor(7)
hThetaEff3P.SetLineWidth(2)
hThetaEff3P.SetTitle('Three-Prong Reconstruction Efficiency vs Theta')
hThetaEff3P.GetXaxis().SetTitle('True Visible #theta [rad]')
hThetaEff3P.GetYaxis().SetTitle('#epsilon')
hThetaEff3P.SetStats(0)
hists.append(hThetaEff3P)

hPhiEff3P = hLinkedMCTauPhi3P.Clone('3p_phi_eff')
hPhiEff3P.Divide(hPhiEff3P, hTauVisPhiTrue3P, 1, 1, 'B')
hPhiEff3P.SetLineColor(418)
hPhiEff3P.SetLineWidth(2)
hPhiEff3P.SetTitle('Three-Prong Reconstruction Efficiency vs Phi')
hPhiEff3P.GetXaxis().SetTitle('True Visible #phi [rad]')
hPhiEff3P.GetYaxis().SetTitle('#epsilon')
hPhiEff3P.SetStats(0)
hists.append(hPhiEff3P)

# Create reco-one-prong efficiency hists
hPtEff1PReco = hLinkedMCTauPt1PReco.Clone('reco_1p_pt_eff')
hPtEff1PReco.Divide(hPtEff1PReco, hTauVisPtTrue1P, 1, 1, 'B')
hPtEff1PReco.SetLineColor(6)
hPtEff1PReco.SetLineWidth(2)
hPtEff1PReco.SetTitle('Reco-One-Prong Reconstruction Efficiency vs Pt')
hPtEff1PReco.GetXaxis().SetTitle('True Visible Pt [GeV/c]')
hPtEff1PReco.GetYaxis().SetTitle('#epsilon')
hPtEff1PReco.SetStats(0)
hists.append(hPtEff1PReco)

hThetaEff1PReco = hLinkedMCTauTheta1PReco.Clone('reco_1p_theta_eff')
hThetaEff1PReco.Divide(hThetaEff1PReco, hTauVisThetaTrue1P, 1, 1, 'B')
hThetaEff1PReco.SetLineColor(7)
hThetaEff1PReco.SetLineWidth(2)
hThetaEff1PReco.SetTitle('Reco-One-Prong Reconstruction Efficiency vs Theta')
hThetaEff1PReco.GetXaxis().SetTitle('True Visible #theta [rad]')
hThetaEff1PReco.GetYaxis().SetTitle('#epsilon')
hThetaEff1PReco.SetStats(0)
hists.append(hThetaEff1PReco)

hPhiEff1PReco = hLinkedMCTauPhi1PReco.Clone('reco_1p_phi_eff')
hPhiEff1PReco.Divide(hPhiEff1PReco, hTauVisPhiTrue1P, 1, 1, 'B')
hPhiEff1PReco.SetLineColor(418)
hPhiEff1PReco.SetLineWidth(2)
hPhiEff1PReco.SetTitle('Reco-One-Prong Reconstruction Efficiency vs Phi')
hPhiEff1PReco.GetXaxis().SetTitle('True Visible #phi [rad]')
hPhiEff1PReco.GetYaxis().SetTitle('#epsilon')
hPhiEff1PReco.SetStats(0)
hists.append(hPhiEff1PReco)

# Create reco-three-prong efficiency hists
hPtEff3PReco = hLinkedMCTauPt3PReco.Clone('reco_3p_pt_eff')
hPtEff3PReco.Divide(hPtEff3PReco, hTauVisPtTrue3P, 1, 1, 'B')
hPtEff3PReco.SetLineColor(6)
hPtEff3PReco.SetLineWidth(2)
hPtEff3PReco.SetTitle('Reco-Three-Prong Reconstruction Efficiency vs Pt')
hPtEff3PReco.GetXaxis().SetTitle('True Visible Pt [GeV/c]')
hPtEff3PReco.GetYaxis().SetTitle('#epsilon')
hPtEff3PReco.SetStats(0)
hists.append(hPtEff3PReco)

hThetaEff3PReco = hLinkedMCTauTheta3PReco.Clone('reco_3p_theta_eff')
hThetaEff3PReco.Divide(hThetaEff3PReco, hTauVisThetaTrue3P, 1, 1, 'B')
hThetaEff3PReco.SetLineColor(7)
hThetaEff3PReco.SetLineWidth(2)
hThetaEff3PReco.SetTitle('Reco-Three-Prong Reconstruction Efficiency vs Theta')
hThetaEff3PReco.GetXaxis().SetTitle('True Visible #theta [rad]')
hThetaEff3PReco.GetYaxis().SetTitle('#epsilon')
hThetaEff3PReco.SetStats(0)
hists.append(hThetaEff3PReco)

hPhiEff3PReco = hLinkedMCTauPhi3PReco.Clone('reco_3p_phi_eff')
hPhiEff3PReco.Divide(hPhiEff3PReco, hTauVisPhiTrue3P, 1, 1, 'B')
hPhiEff3PReco.SetLineColor(418)
hPhiEff3PReco.SetLineWidth(2)
hPhiEff3PReco.SetTitle('Reco-Three-Prong Reconstruction Efficiency vs Phi')
hPhiEff3PReco.GetXaxis().SetTitle('True Visible #phi [rad]')
hPhiEff3PReco.GetYaxis().SetTitle('#epsilon')
hPhiEff3PReco.SetStats(0)
hists.append(hPhiEff3PReco)

# Create fake one-prong hists
hPtFake1P = hPtEff1P.Clone('fake_1p_pt')
hPtFake1P.Add(hPtEff1PReco, -1)
hPtFake1P.SetLineColor(6)
hPtFake1P.SetLineWidth(2)
hPtFake1P.SetTitle('Fake One-Prong Rate vs Pt')
hPtFake1P.GetXaxis().SetTitle('True Visible Pt [GeV/c]')
hPtFake1P.GetYaxis().SetTitle('#epsilon')
hPtFake1P.SetStats(0)
hists.append(hPtFake1P)

hThetaFake1P = hThetaEff1P.Clone('fake_1p_theta')
hThetaFake1P.Add(hThetaEff1PReco, -1)
hThetaFake1P.SetLineColor(7)
hThetaFake1P.SetLineWidth(2)
hThetaFake1P.SetTitle('Fake One-Prong Rate vs Theta')
hThetaFake1P.GetXaxis().SetTitle('True Visible #theta [rad]')
hThetaFake1P.GetYaxis().SetTitle('#epsilon')
hThetaFake1P.SetStats(0)
hists.append(hThetaFake1P)

hPhiFake1P = hPhiEff1P.Clone('fake_1p_phi')
hPhiFake1P.Add(hPhiEff1PReco, -1)
hPhiFake1P.SetLineColor(418)
hPhiFake1P.SetLineWidth(2)
hPhiFake1P.SetTitle('Fake One-Prong Rate vs Phi')
hPhiFake1P.GetXaxis().SetTitle('True Visible #phi [rad]')
hPhiFake1P.GetYaxis().SetTitle('#epsilon')
hPhiFake1P.SetStats(0)
hists.append(hPhiFake1P)

# Create fake three-prong hists
hPtFake3P = hPtEff3P.Clone('fake_3p_pt')
hPtFake3P.Add(hPtEff3PReco, -1)
hPtFake3P.SetLineColor(6)
hPtFake3P.SetLineWidth(2)
hPtFake3P.SetTitle('Fake Three-Prong Rate vs Pt')
hPtFake3P.GetXaxis().SetTitle('True Visible Pt [GeV/c]')
hPtFake3P.GetYaxis().SetTitle('#epsilon')
hPtFake3P.SetStats(0)
hists.append(hPtFake3P)

hThetaFake3P = hThetaEff3P.Clone('fake_3p_theta')
hThetaFake3P.Add(hThetaEff3PReco, -1)
hThetaFake3P.SetLineColor(7)
hThetaFake3P.SetLineWidth(2)
hThetaFake3P.SetTitle('Fake Three-Prong Rate vs Theta')
hThetaFake3P.GetXaxis().SetTitle('True Visible #theta [rad]')
hThetaFake3P.GetYaxis().SetTitle('#epsilon')
hThetaFake3P.SetStats(0)
hists.append(hThetaFake3P)

hPhiFake3P = hPhiEff3P.Clone('fake_3p_phi')
hPhiFake3P.Add(hPhiEff3PReco, -1)
hPhiFake3P.SetLineColor(418)
hPhiFake3P.SetLineWidth(2)
hPhiFake3P.SetTitle('Fake Three-Prong Rate vs Phi')
hPhiFake3P.GetXaxis().SetTitle('True Visible #phi [rad]')
hPhiFake3P.GetYaxis().SetTitle('#epsilon')
hPhiFake3P.SetStats(0)
hists.append(hPhiFake3P)

# Create 1P Pion Eff Hists
hPiPtEff1P = hPiPtMatched1P.Clone('1p_pi_pt_eff')
hPiPtEff1P.Divide(hPiPtEff1P, hPiPtTrue1P, 1, 1, 'B')
hPiPtEff1P.SetLineColor(6)
hPiPtEff1P.SetLineWidth(2)
hPiPtEff1P.SetTitle('One-Prong Pion Reconstruction Efficiency vs Pt')
hPiPtEff1P.GetXaxis().SetTitle('True Visible Pt [GeV/c]')
hPiPtEff1P.GetYaxis().SetTitle('#epsilon')
hPiPtEff1P.SetStats(0)
hists.append(hPiPtEff1P)

hPiThetaEff1P = hPiThetaMatched1P.Clone('1p_pi_theta_eff')
hPiThetaEff1P.Divide(hPiThetaEff1P, hPiThetaTrue1P, 1, 1, 'B')
hPiThetaEff1P.SetLineColor(7)
hPiThetaEff1P.SetLineWidth(2)
hPiThetaEff1P.SetTitle('One-Prong Pion Reconstruction Efficiency vs Theta')
hPiThetaEff1P.GetXaxis().SetTitle('True Visible #theta [rad]')
hPiThetaEff1P.GetYaxis().SetTitle('#epsilon')
hPiThetaEff1P.SetStats(0)
hists.append(hPiThetaEff1P)

hPiPhiEff1P = hPiPhiMatched1P.Clone('1p_pi_phi_eff')
hPiPhiEff1P.Divide(hPiPhiEff1P, hPiPhiTrue1P, 1, 1, 'B')
hPiPhiEff1P.SetLineColor(418)
hPiPhiEff1P.SetLineWidth(2)
hPiPhiEff1P.SetTitle('One-Prong Pion Reconstruction Efficiency vs Phi')
hPiPhiEff1P.GetXaxis().SetTitle('True Visible #phi [rad]')
hPiPhiEff1P.GetYaxis().SetTitle('#epsilon')
hPiPhiEff1P.SetStats(0)
hists.append(hPiPhiEff1P)

# Create 3P Pion Eff Hists
hPiPtEff3P = hPiPtMatched3P.Clone('3p_pi_pt_eff')
hPiPtEff3P.Divide(hPiPtEff3P, hPiPtTrue3P, 1, 1, 'B')
hPiPtEff3P.SetLineColor(6)
hPiPtEff3P.SetLineWidth(2)
hPiPtEff3P.SetTitle('Three-Prong Pion Reconstruction Efficiency vs Pt')
hPiPtEff3P.GetXaxis().SetTitle('True Visible Pt [GeV/c]')
hPiPtEff3P.GetYaxis().SetTitle('#epsilon')
hPiPtEff3P.SetStats(0)
hists.append(hPiPtEff3P)

hPiThetaEff3P = hPiThetaMatched3P.Clone('3p_pi_theta_eff')
hPiThetaEff3P.Divide(hPiThetaEff3P, hPiThetaTrue3P, 1, 1, 'B')
hPiThetaEff3P.SetLineColor(7)
hPiThetaEff3P.SetLineWidth(2)
hPiThetaEff3P.SetTitle('Three-Prong Pion Reconstruction Efficiency vs Theta')
hPiThetaEff3P.GetXaxis().SetTitle('True Visible #theta [rad]')
hPiThetaEff3P.GetYaxis().SetTitle('#epsilon')
hPiThetaEff3P.SetStats(0)
hists.append(hPiThetaEff3P)

hPiPhiEff3P = hPiPhiMatched3P.Clone('3p_pi_phi_eff')
hPiPhiEff3P.Divide(hPiPhiEff3P, hPiPhiTrue3P, 1, 1, 'B')
hPiPhiEff3P.SetLineColor(438)
hPiPhiEff3P.SetLineWidth(2)
hPiPhiEff3P.SetTitle('Three-Prong Pion Reconstruction Efficiency vs Phi')
hPiPhiEff3P.GetXaxis().SetTitle('True Visible #phi [rad]')
hPiPhiEff3P.GetYaxis().SetTitle('#epsilon')
hPiPhiEff3P.SetStats(0)
hists.append(hPiPhiEff3P)

# Write to output file
output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
    hist.Write()
output_file.Close()

# Draw hists and save as PNG
for hist in hists:
    filename = hist.GetTitle() + '.png'
    canvas = TCanvas()
    hist.Draw()
    canvas.SaveAs(filename)

with open('three_prong_types.txt', 'w') as file:
    file.write('Three-Prong Daughter Types')
    for key, value, in three_prong_types.items():
        print(f'Type: {key}, Total Number: {value}\n', file=file)
