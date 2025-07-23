from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt

from tau_mc_link import getLinkedMCTau, getVisibleProperties, getDecayMode, getNRecoQPis

def getEta(mom):

    p = math.sqrt(mom[0]**2 + mom[1]**2 + mom[2]**2)
    pz = mom[2]
    theta = math.acos(pz/p)
    return -math.log(math.tan(theta/2))
    
def getPhi(mom):
    return math.atan2(mom[1], mom[0])

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='reco_output.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='pi_ana.root')

args = parser.parse_args()

# Initialize histograms
hists = []

hPiMaxConeTrue = TH1F('true_pi_max_angle', 'True Maximum Search Cone Angle', 50, 0, 0.25)
hPiMaxConeTrue.GetXaxis().SetTitle('#theta_{max} [rad]')
hists.append(hPiMaxConeTrue)

hPiMaxConeReco = TH1F('reco_pi_max_angle', 'Reconstructed Maximum Search Cone Angle', 50, 0, 0.25)
hPiMaxConeReco.GetXaxis().SetTitle('#theta_{max} [rad]')
hists.append(hPiMaxConeReco)

hPiPtTrue1P0N = TH1F('true_pi_pt_1p0n', 'True 1P0N Pion Pt', 100, 0, 320)
hists.append(hPiPtTrue1P0N)

hPiPtReco1P0N = TH1F('reco_pi_pt_1p0n', 'Reco 1P0N Pion Pt', 100, 0, 2000)
hists.append(hPiPtReco1P0N)

hPiPtTrue3P0N = TH1F('true_pi_pt_3p0n', 'True 3P0N Pion Pt', 100, 0, 320)
hists.append(hPiPtTrue3P0N)

hPiPtReco3P0N = TH1F('reco_pi_pt_3p0n', 'Reco 3P0N Pion Pt', 100, 0, 500)
hists.append(hPiPtReco3P0N)

hPiEtaRatio1P0N = TH1F('pi_eta_ratio_1p0n', '1P0N Pion Eta Ratio', 100, 0, 2.5)
hists.append(hPiEtaRatio1P0N)

hPiPhiRatio1P0N = TH1F('pi_phi_ratio_1p0n', '1P0N Pion Phi Ratio', 100, 0, 2)
hists.append(hPiPhiRatio1P0N)

hNTruePis1P0N = TH1F('n_true_pis_1p0n', 'Number of True Pions (1P0N)', 65, 0, 65)
hists.append(hNTruePis1P0N)

hNTruePis3P0N = TH1F('n_true_pis_3p0n', 'Number of True Pions (3P0N)', 65, 0, 65)
hists.append(hNTruePis3P0N)

# Keep track of max cone angles
pi_max_cone_reco = []
pi_max_cone_true = []

# Keep track of seed pt
pi_seed_pt_reco = []
pi_seed_pt_true = []

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

n_big_pt = 0
    
# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        # Get collections
        pfos = event.getCollection('PandoraPFOs')
        mcParticles = event.getCollection('MCParticle')

        # Loop through MC particles
        for mcParticle in mcParticles:

            # Tag MC taus
            pdg = abs(mcParticle.getPDG())
            if pdg == 15:

                # Get tau decay mode
                decayMode = getDecayMode(mcParticle)

                if decayMode == 0:
                    daughters = mcParticle.getDaughters()
                    energy_max = 0
                    for daughter in daughters:
                        if (abs(daughter.getPDG()) == 211):
                            mom = daughter.getMomentum()
                            pt = math.sqrt(mom[0]**2 + mom[1]**2)
                            hPiPtTrue1P0N.Fill(pt)
                            energy = daughter.getEnergy()
                            if energy > energy_max:
                                energy_max = energy
                                eta_true = getEta(mom)
                                phi_true = getPhi(mom)

                    energy_max = 0
                    for pfo in pfos:
                        if (abs(pfo.getType()) == 211):
                            mom = pfo.getMomentum()
                            pt = math.sqrt(mom[0]**2 + mom[1]**2)
                            hPiPtReco1P0N.Fill(pt)
                            energy = pfo.getEnergy()
                            if energy > energy_max:
                                energy_max = energy
                                eta_reco = getEta(mom)
                                phi_reco = getPhi(mom)

                    eta_ratio = eta_reco/eta_true
                    phi_ratio = phi_reco/phi_true

                    hPiEtaRatio1P0N.Fill(eta_ratio)
                    hPiPhiRatio1P0N.Fill(phi_ratio)

                    n_pis_true_1p0n = 0
                    for mc_particle in mcParticles:
                        if (abs(mc_particle.getPDG()) == 211):
                            n_pis_true_1p0n += 1
                    hNTruePis1P0N.Fill(n_pis_true_1p0n)
                    

                if decayMode == 4:
                    daughters = mcParticle.getDaughters()
                    for daughter in daughters:
                        if (abs(daughter.getPDG()) == 211):
                            mom = daughter.getMomentum()
                            pt = math.sqrt(mom[0]**2 + mom[1]**2)
                            hPiPtTrue3P0N.Fill(pt)

                    for pfo in pfos:
                        if (abs(pfo.getType()) == 211):
                            mom = pfo.getMomentum()
                            pt = math.sqrt(mom[0]**2 + mom[1]**2)
                            hPiPtReco3P0N.Fill(pt)
                            if pt > 500:
                                n_big_pt += 1

                    n_pis_true_3p0n = 0
                    for mc_particle in mcParticles:
                        if (abs(mc_particle.getPDG()) == 211):
                            n_pis_true_3p0n += 1
                    hNTruePis3P0N.Fill(n_pis_true_3p0n)

                if decayMode == 4:

                    pis = []
                    for pfo in pfos:
                        if abs(pfo.getType()) == 211:
                            pis.append(pfo)

                    if len(pis) > 1:
                        # Tag seed
                        pt_seed = 0
                        pi_seed = None
                        for pi in pis:
                            pt = math.sqrt(pi.getMomentum()[0]**2+pi.getMomentum()[1]**2)
                            if pt > pt_seed:
                                pt_seed = pt
                                pi_seed = pi

                        pi_seed_pt_reco.append(pt_seed)

                        px_seed = pi_seed.getMomentum()[0]
                        py_seed = pi_seed.getMomentum()[1]
                        pz_seed = pi_seed.getMomentum()[2]
                        pis.remove(pi_seed)
                        
                        max_angle_reco = 0
                        for pi in pis:
                            px_pi = pi.getMomentum()[0]
                            py_pi = pi.getMomentum()[1]
                            pz_pi = pi.getMomentum()[2]
                            angle = math.acos((px_seed*px_pi+py_seed*py_pi+pz_seed*pz_pi)/(math.sqrt(px_seed**2+py_seed**2+pz_seed**2)*math.sqrt(px_pi**2+py_pi**2+pz_pi**2)))
                            if angle > max_angle_reco:
                                max_angle_reco = angle

                        hPiMaxConeReco.Fill(max_angle_reco)
                        pi_max_cone_reco.append(max_angle_reco)
                            
                    pis_true = []
                    daughters = mcParticle.getDaughters()
                    for daughter in daughters:
                        if abs(daughter.getPDG()) == 211:
                            pis_true.append(daughter)
                            
                    if len(pis_true) > 1:
                        # Tag seed
                        pt_seed = 0
                        pi_seed = None
                        for pi_true in pis_true:
                            pt = math.sqrt(pi_true.getMomentum()[0]**2+pi_true.getMomentum()[1]**2)
                            if pt > pt_seed:
                                pt_seed = pt
                                pi_seed = pi_true

                        pi_seed_pt_true.append(pt_seed)

                        px_seed = pi_seed.getMomentum()[0]
                        py_seed = pi_seed.getMomentum()[1]
                        pz_seed = pi_seed.getMomentum()[2]
                        pis_true.remove(pi_seed)
                        
                        max_angle_true = 0
                        for pi_true in pis_true:
                            px_pi = pi_true.getMomentum()[0]
                            py_pi = pi_true.getMomentum()[1]
                            pz_pi = pi_true.getMomentum()[2]
                            angle = math.acos((px_seed*px_pi+py_seed*py_pi+pz_seed*pz_pi)/(math.sqrt(px_seed**2+py_seed**2+pz_seed**2)*math.sqrt(px_pi**2+py_pi**2+pz_pi**2)))
                            if angle > max_angle_true:
                                max_angle_true = angle

                        hPiMaxConeTrue.Fill(max_angle_true)
                        pi_max_cone_true.append(max_angle_true)

    # Close file
    reader.close()

print(f'Number of 3P0N pions with pt > 500: {n_big_pt}')

plt.figure()

plt.clf()
plt.scatter(pi_seed_pt_reco, pi_max_cone_reco, c='blue')
plt.xlim(0, 250)
plt.xlabel('Reco Seed Pt [GeV/c]')
plt.ylabel('Reco Max Angle [rad]')
plt.savefig('seed_cone_reco_250.png')

plt.clf()
plt.scatter(pi_seed_pt_reco, pi_max_cone_reco, c='blue')
plt.xlim(0, 2000)
plt.xlabel('Reco Seed Pt [GeV/c]')
plt.ylabel('Reco Max Angle [rad]')
plt.savefig('seed_cone_reco_2000.png')

plt.clf()
plt.scatter(pi_seed_pt_reco, pi_max_cone_reco, c='blue')
plt.xlabel('Reco Seed Pt [GeV/c]')
plt.ylabel('Reco Max Angle [rad]')
plt.savefig('seed_cone_reco_max.png')

plt.clf()
plt.scatter(pi_seed_pt_true, pi_max_cone_true, c='blue')
plt.xlabel('True Seed Pt [GeV/c]')
plt.ylabel('True Max Angle [rad]')
plt.savefig('seed_cone_true.png')
    
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
