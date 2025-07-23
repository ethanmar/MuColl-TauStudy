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

# Script to analyze individual tau decay modes

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

args = parser.parse_args()

# Check if input file is a directory or a single file    
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

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        taus_passed = event.getCollection('RecoTaus')
        taus_nparticles = event.getCollection('RecoTausNParticles')
        taus_nchargedtrks = event.getCollection('RecoTausNChargedTrks')
        taus_merge = event.getCollection('RecoTausMerge')
        taus_invmass = event.getCollection('RecoTausInvMass')
        taus_isoenergy = event.getCollection('RecoTausIsoEnergy')
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')
        tauRecoLink = event.getCollection('TauPFOLink')
        recoMCLink = event.getCollection('RecoMCTruthLink')

        # Instantiate relation navigators to parse tauReco and RecoMC links
        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        # Merge all failed tau collections into one list
        taus_failed = taus_nparticles + taus_nchargedtrks + taus_merge + taus_invmass + taus_isoenergy

        for tau_passed in taus_passed:
            true_tau = getLinkedMCTau(tau_passed, relationNavigatorTau, relationNavigatorRecoMC)
            decay_mode = getDecayMode(true_tau)
            if decay_mode == 4:
                n_unique_pis, n_duplicate_pis = getNRecoPis(tau_passed, relationNavigatorRecoMC)
                if (n_unique_pis == 3 and n_duplicate_pis == 0):
                    vis_props = getVisibleProperties(true_tau)
                    vis_pt = getPt(vis_props)
                    vis_theta = getTheta(vis_props)
                    vis_phi = getPhi(vis_props)
                    if (vis_theta < 0.577 or vis_theta > 2.56):
                        if ((65 <= vis_pt and vis_pt <= 128) or (193 <= vis_pt and vis_pt <= 224)):
                            print(f'\nPassed Tau')
                            print(f'Visible Theta: {vis_theta}, Visible Pt: {vis_pt}')
                            pfos = tau_passed.getParticles()
                            for pfo in pfos:
                                if (abs(pfo.getType()) == 211):
                                    mom = pfo.getMomentum()
                                    pt = math.sqrt(mom[1]**2 + mom[0]**2)
                                    print(f'Pion Reco Pt: {pt}')
        for tau_failed in taus_failed:
            true_tau = getLinkedMCTau(tau_failed, relationNavigatorTau, relationNavigatorRecoMC)
            decay_mode = getDecayMode(true_tau)
            if decay_mode == 4:
                n_unique_pis, n_duplicate_pis = getNRecoPis(tau_failed, relationNavigatorRecoMC)
                if (n_unique_pis == 3 and n_duplicate_pis == 0):
                    vis_props = getVisibleProperties(true_tau)
                    vis_pt = getPt(vis_props)
                    vis_theta = getTheta(vis_props)
                    vis_phi = getPhi(vis_props)
                    if (vis_theta < 0.577 or vis_theta > 2.56):
                        if ((65 <= vis_pt and vis_pt <= 128) or (193 <= vis_pt and vis_pt <= 224)):
                            print(f'\nFailed Tau')
                            print(f'Visible Theta: {vis_theta}, Visible Pt: {vis_pt}')
                            pfos = tau_failed.getParticles()
                            for pfo in pfos:
                                if (abs(pfo.getType()) == 211):
                                    mom = pfo.getMomentum()
                                    pt = math.sqrt(mom[1]**2 + mom[0]**2)
                                    print(f'Pion Reco Pt: {pt}')
    reader.close()
