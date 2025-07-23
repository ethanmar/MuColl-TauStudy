from pyLCIO import IOIMPL, EVENT, UTIL
from argparse import ArgumentParser
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

from tau_mc_link import getDecayMode, getLinkedMCTau

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

# Keep track of true tau label
true_tau_label = []

# Keep track of reco tau label
reco_tau_label = []

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
        taus_invmass = event.getCollection('RecoTausInvMass')
        taus_isoenergy = event.getCollection('RecoTausIsoEnergy')
        pfos = event.getCollection('PandoraPFOs')
        tauRecoLink = event.getCollection('TauPFOLink')
        recoMCLink = event.getCollection('RecoMCTruthLink')
        
        # Instantiate relation navigators to parse tauReco and RecoMC links
        relationNavigatorTau = UTIL.LCRelationNavigator(tauRecoLink)
        relationNavigatorRecoMC = UTIL.LCRelationNavigator(recoMCLink)
        
        # Reco taus that passed
        for tau in taus:

            # Get true label
            true_tau = getLinkedMCTau(tau, relationNavigatorTau, relationNavigatorRecoMC)
            true_decay_mode = getDecayMode(true_tau)
            if (true_decay_mode == 0):
                true_tau_label.append(1)
            elif (true_decay_mode == 4):
                true_tau_label.append(3)

            # Get reco label
            if (decay_mode == 0 or decay_mode == 4):
                pfos = tau.getParticles()
                n_pis = 0
                for pfo in pfos:
                    if (abs(pfo.getType()) == 211):
                        n_pis += 1
                reco_tau_label.append(n_pis)
                
        '''
        # Reco taus that failed number of charged tracks
        for tau_nchargedtrks in taus_nchargedtrks:

            # Get true label
            true_tau = getLinkedMCTau(tau_nchargedtrks, relationNavigatorTau, relationNavigatorRecoMC)
            true_decay_mode = getDecayMode(true_tau)
            if (true_decay_mode != 0 and true_decay_mode != 4):
                continue
            elif (true_decay_mode == 0):
                true_tau_label.append(1)
            elif (true_decay_mode == 4):
                true_tau_label.append(3)

            # Get reco label
            pfos = tau_nchargedtrks.getParticles()
            n_pis = 0
            for pfo in pfos:
                if (abs(pfo.getType()) == 211):
                    n_pis += 1
            reco_tau_label.append(n_pis)
            
        # Reco taus that failed invariant mass
        for tau_invmass in taus_invmass:

            # Get true label
            true_tau = getLinkedMCTau(tau_invmass, relationNavigatorTau, relationNavigatorRecoMC)
            true_decay_mode = getDecayMode(true_tau)
            if (true_decay_mode != 0 and true_decay_mode != 4):
                continue
            elif (true_decay_mode == 0):
                true_tau_label.append(1)
            elif (true_decay_mode == 4):
                true_tau_label.append(3)

            # Get reco label
            pfos = tau_invmass.getParticles()
            n_pis = 0
            for pfo in pfos:
                if (abs(pfo.getType()) == 211):
                    n_pis += 1
            reco_tau_label.append(n_pis)

        # Reco taus that failed merge
        for tau_merge in taus_merge:

            # Get true label
            true_tau = getLinkedMCTau(tau_merge, relationNavigatorTau, relationNavigatorRecoMC)
            true_decay_mode = getDecayMode(true_tau)
            if (true_decay_mode != 0 and true_decay_mode != 4):
                continue
            elif (true_decay_mode == 0):
                true_tau_label.append(1)
            elif (true_decay_mode == 4):
                true_tau_label.append(3)

            # Get reco label
            pfos = tau_merge.getParticles()
            n_pis = 0
            for pfo in pfos:
                if (abs(pfo.getType()) == 211):
                    n_pis += 1
            reco_tau_label.append(n_pis)

        # Reco taus that failed number of particles
        for tau_nparticles in taus_nparticles:

            # Get true label
            true_tau = getLinkedMCTau(tau_nparticles, relationNavigatorTau, relationNavigatorRecoMC)
            true_decay_mode = getDecayMode(true_tau)
            if (true_decay_mode != 0 and true_decay_mode != 4):
                continue
            elif (true_decay_mode == 0):
                true_tau_label.append(1)
            elif (true_decay_mode == 4):
                true_tau_label.append(3)

            # Get reco label
            pfos = tau_nparticles.getParticles()
            n_pis = 0
            for pfo in pfos:
                if (abs(pfo.getType()) == 211):
                    n_pis += 1
            reco_tau_label.append(n_pis)

        # Reco taus that failed isolation energy
        for tau_isoenergy in taus_isoenergy:

            # Get true label
            true_tau = getLinkedMCTau(tau_isoenergy, relationNavigatorTau, relationNavigatorRecoMC)
            true_decay_mode = getDecayMode(true_tau)
            if (true_decay_mode != 0 and true_decay_mode != 4):
                continue
            elif (true_decay_mode == 0):
                true_tau_label.append(1)
            elif (true_decay_mode == 4):
                true_tau_label.append(3)

            # Get reco label
            pfos = tau_isoenergy.getParticles()
            n_pis = 0
            for pfo in pfos:
                if (abs(pfo.getType()) == 211):
                    n_pis += 1
            reco_tau_label.append(n_pis)
        '''

    # Close file
    reader.close()

# Create confusion matrix
# labels=['1-Prong', '2-Prong', '3-Prong', '4-Prong']
disp = ConfusionMatrixDisplay.from_predictions(reco_tau_label, true_tau_label, normalize='pred', include_values=True, cmap='Blues', colorbar=True)
disp.ax_.set_ylabel(r'Number of Reconstructed $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
disp.ax_.set_xlabel(r'Number of True $\pi^\pm$s (0 $\pi^0$s)', fontsize=12)
disp.ax_.invert_yaxis()
plt.tight_layout()
plt.savefig('tau_id_confusion.png')
plt.show()

