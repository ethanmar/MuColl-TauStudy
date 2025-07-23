from pyLCIO import IOIMPL, EVENT, UTIL
from argparse import ArgumentParser
import os
import numpy as np
import math
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

from tau_mc_link import getDecayMode, getLinkedMCTau

def getEta(mom):

    p = math.sqrt(mom[0]**2 + mom[1]**2 + mom[2]**2)
    pz = mom[2]
    theta = math.acos(pz/p)
    return -math.log(math.tan(theta/2))
    
def getPhi(mom):
    return math.atan2(mom[1], mom[0])

def getDeltaR(mom_1, mom_2):
    eta_1 = getEta(mom_1)
    phi_1 = getPhi(mom_1)
    
    eta_2 = getEta(mom_2)
    phi_2 = getPhi(mom_2)

    delta_eta = eta_1 - eta_2
    delta_phi = phi_1 - phi_2

    if (delta_phi > math.pi):
        delta_phi -= 2*math.pi
    elif (delta_phi < math.pi):
        delta_phi += 2*math.pi

    return math.sqrt(delta_eta**2 + delta_phi**2)

'''
def getRecoPDG(mcParticle, pfos):

    # Get mcParticle momentum
    mom_mc = mcParticle.getMomentum()

    # Loop through pfos
    energy_max = 0
    pdg_reco = 0
    for pfo in pfos:
        
        # Get pfo momentum
        mom_pfo = pfo.getMomentum()

        # Calculate deltaR
        deltaR = getDeltaR(mom_mc, mom_pfo)
        if (abs(pfo.getType()) == 211):
            print(f'DeltaR of Pion: {deltaR}')

        # Check within deltaR cone of 0.1
        if (deltaR < 0.1):
            energy = pfo.getEnergy()
            if energy > energy_max:
                energy_max = energy
                pdg_reco = abs(pfo.getType())

    return pdg_reco
'''

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

    '''
# Keep track of true pdg
true_pdg_1P0N = []
true_pdg_3P0N = []

# Keep track of reco pdg
reco_pdg_1P0N = []
reco_pdg_3P0N = []
'''

# Keep track of deltaR
deltaR_1p0n = []
deltaR_3p0n = []

# Keep track of etas
eta_true_1p0n = []
eta_reco_1p0n = []
eta_true_3p0n = []
eta_reco_3p0n = []

# Keep track of phis
phi_true_1p0n = []
phi_reco_1p0n = []
phi_true_3p0n = []
phi_reco_3p0n = []


# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')
        recoMCLink = event.getCollection('RecoMCTruthLink')
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)

        for mcParticle in mcParticles:
            if abs(mcParticle.getPDG()) == 15:
                decay_mode = getDecayMode(mcParticle)

        if decay_mode == 0 or decay_mode == 4:
            for mcParticle in mcParticles:
                true_pdg = abs(mcParticle.getPDG())
                if true_pdg != 211:
                    continue
                true_mom = mcParticle.getMomentum()
                true_phi = getPhi(true_mom)
                true_eta = getEta(true_mom)
                for pfo in pfos:
                    reco_pdg = abs(pfo.getType())
                    if reco_pdg == 211:
                        reco_mom = pfo.getMomentum()
                        reco_phi = getPhi(reco_mom)
                        reco_eta = getEta(reco_mom)
                        deltaR = getDeltaR(true_mom, reco_mom)
                        if decay_mode == 0:
                            deltaR_1p0n.append(deltaR)
                            phi_true_1p0n.append(true_phi)
                            phi_reco_1p0n.append(reco_phi)
                            eta_true_1p0n.append(true_eta)
                            eta_reco_1p0n.append(reco_eta)
                        elif decay_mode == 4:
                            deltaR_3p0n.append(deltaR)
                            phi_true_3p0n.append(true_phi)
                            phi_reco_3p0n.append(reco_phi)
                            eta_true_3p0n.append(true_eta)
                            eta_reco_3p0n.append(reco_eta)
                '''
                reco_pdg = getRecoPDG(mcParticle, pfos)
                if decay_mode == 0:
                    reco_pdg_1P0N.append(reco_pdg)
                    true_pdg_1P0N.append(true_pdg)
                elif decay_mode == 4:
                    reco_pdg_3P0N.append(reco_pdg)
                    true_pdg_3P0N.append(true_pdg)     
                '''
                
    reader.close()
                    
plt.figure(figsize=(10,10))

plt.clf()
plt.hist(deltaR_1p0n, bins=25, color='blue')
plt.xlabel('DeltaR (1P0N)')
plt.savefig('deltaR_1p0n.png')

plt.clf()
plt.hist(deltaR_3p0n, bins=25, color='blue')
plt.xlabel('DeltaR (3P0N)')
plt.savefig('deltaR_3p0n.png')

plt.clf()
plt.hist(phi_true_1p0n, bins=25, color='blue')
plt.xlabel('True Phi (1P0N)')
plt.savefig('phi_true_1p0n.png')

plt.clf()
plt.hist(phi_reco_1p0n, bins=25, color='blue')
plt.xlabel('Reco Phi (1P0N)')
plt.savefig('phi_reco_1p0n.png')

plt.clf()
plt.hist(phi_true_3p0n, bins=25, color='blue')
plt.xlabel('True Phi (3P0N)')
plt.savefig('phi_true_3p0n.png')

plt.clf()
plt.hist(phi_reco_3p0n, bins=25, color='blue')
plt.xlabel('Reco Phi (3P0N)')
plt.savefig('phi_reco_3p0n.png')

plt.clf()
plt.hist(eta_true_1p0n, bins=25, color='blue')
plt.xlabel('True Eta (1P0N)')
plt.savefig('eta_true_1p0n.png')

plt.clf()
plt.hist(eta_reco_1p0n, bins=25, color='blue')
plt.xlabel('Reco Eta (1P0N)')
plt.savefig('eta_reco_1p0n.png')

plt.clf()
plt.hist(eta_true_3p0n, bins=25, color='blue')
plt.xlabel('True Eta (3P0N)')
plt.savefig('eta_true_3p0n.png')

plt.clf()
plt.hist(eta_reco_3p0n, bins=25, color='blue')
plt.xlabel('Reco Eta (3P0N)')
plt.savefig('eta_reco_3p0n.png')

'''
# Create 1P0N confusion matrix
plt.clf()
disp_1P0N = ConfusionMatrixDisplay.from_predictions(reco_pdg_1P0N, true_pdg_1P0N, normalize='pred', include_values=True, cmap='Blues', colorbar=True)
disp_1P0N.ax_.set_ylabel('Reconstructed Particle Type (1P0N)', fontsize=12)
disp_1P0N.ax_.set_xlabel('True Particle Type (1P0N)', fontsize=12)
disp_1P0N.ax_.invert_yaxis()
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('pdg_id_1P0N_confusion.png')

# Create 3P0N confusion matrix
plt.clf()
disp_3P0N = ConfusionMatrixDisplay.from_predictions(reco_pdg_3P0N, true_pdg_3P0N, normalize='pred', include_values=True, cmap='Blues', colorbar=True)
disp_3P0N.ax_.set_ylabel('Reconstructed Particle Type (3P0N)', fontsize=12)
disp_3P0N.ax_.set_xlabel('True Particle Type (3P0N)', fontsize=12)
disp_3P0N.ax_.invert_yaxis()
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('pdg_id_3P0N_confusion.png')
'''
