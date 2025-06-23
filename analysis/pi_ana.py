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

    # Close file
    reader.close()
    
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
