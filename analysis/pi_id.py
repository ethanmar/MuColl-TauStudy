from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
import os
import numpy as np

def getMatchedMCParticle(pfo, mc_particles):

    # Get pfo observables
    pfo_mom = pfo.getMomentum()
    pfo_px = pfo_mom[0]
    pfo_py = pfo_mom[1]
    pfo_pz = pfo_mom[2]
    pfo_pt = math.sqrt(pfo_px**2 + pfo_py**2)
    pfo_theta = math.atan2(pfo_pt, pfo_pz)
    pfo_phi = math.atan2(pfo_py, pfo_px)
    pfo_energy = pfo.getEnergy()

    # Loop over mc particles
    deltaR_mc_particles = []
    for mc_particle in mc_particles:

        # Get mc particle observables
        mc_mom = mc_particle.getMomentum()
        mc_px = mc_mom[0]
        mc_py = mc_mom[1]
        mc_pz = mc_mom[2]
        mc_pt = math.sqrt(mc_px**2 + mc_py**2)
        mc_theta = math.atan2(mc_pt, mc_pz)
        mc_phi = math.atan2(mc_py, mc_px)

        # Compute deltaR
        deltaTheta = pfo_theta - mc_theta
        deltaPhi = pfo_phi - mc_phi
        while (deltaPhi > math.pi):
            deltaPhi -= 2*math.pi
        while (deltaPhi < -math.pi):
            deltaPhi += 2*math.pi
        deltaR = math.sqrt(deltaTheta**2 + deltaPhi**2)

        # Save deltaR if less than 0.1
        if (deltaR < 0.1):
            deltaR_mc_particles.append(mc_particle)

    # Loop over mc particles with deltaR less than 0.1
    min_energy_diff = 1e6
    matched_mc_particle = None
    for deltaR_mc_particle in deltaR_mc_particles:

        # Compute energy diference
        mc_energy = deltaR_mc_particle.getEnergy()
        energy_diff = abs(pfo_energy - mc_energy)

        if (energy_diff < min_energy_diff):
            min_energy_diff = energy_diff
            matched_mc_particle = deltaR_mc_particle

    return matched_mc_particle

def fromTau(mc_particle):

    from_tau = False
    
    # Loop through parents
    parents = mc_particle.getParents()
    for parent in parents:
        if (abs(parent.getPDG()) == 15):
            from_tau = True

    return from_tau

def getTheta(pfo):
    mom = pfo.getMomentum()
    px = mom[0]
    py = mom[1]
    pz = mom[2]
    pt = math.sqrt(px**2 + py**2)
    theta = math.atan2(pt, pz)
    return theta

def getLayer(r):
    # 50 layers in ECAL barrel (5.35 mm)
    # 50 layers in ECAL endcap (36.3 mm)
    # 75 layers in HCAL barrel (26.5 mm)
    # 75 layers in HCAL endcap (26.5 mm)

# Command line arguments
parser = ArgumentParser()

# Input file(s)
parser.add_argument('--inputFile', type=str, default='reco_output.slcio')

# Output file
parser.add_argument('--outputFile', type=str, default='pi_id.root')

args = parser.parse_args()

# Initialize histograms
hists = []

# Central Barrel

# Centroid R
hPiCentroidRCentralBarrel = TH1F('pi_centroid_r_cent_barrel', 'Pion Centroid Radius (Central Barrel)', 500, 1800, 4200)
hists.append(hPiCentroidRCentralBarrel)

hElecCentroidRCentralBarrel = TH1F('elec_centroid_r_cent_barrel', 'Electron Centroid Radius (Central Barrel)', 500, 1800, 4200)
hists.append(hElecCentroidRCentralBarrel)

# E/(E+H)
hPiEOverEPlusHCentralBarrel = TH1F('pi_e_over_e_plus_h_cent_barrel', 'Pion E/(E+H) (Central Barrel)', 50, 0, 1.1)
hists.append(hPiEOverEPlusHCentralBarrel)

hElecEOverEPlusHCentralBarrel = TH1F('elec_e_over_e_plus_h_cent_barrel', 'Electron E/(E+H) (Central Barrel)', 50, 0, 1.1)
hists.append(hElecEOverEPlusHCentralBarrel)

# Shower start layer
hPiShowerStartLayerCentralBarrel = TH1F('pi_shower_start_layer_cent_barrel', 'Pion Shower Start Layer (Central Barrel)', 21, 0, 20)
hists.append(hPiShowerStartLayerCentralBarrel)

hElecShowerStartLayerCentralBarrel = TH1F('elec_shower_start_layer_cent_barrel', 'Electron Shower Start Layer (Central Barrel)', 21, 0, 20)
hists.append(hElecShowerStartLayerCentralBarrel)

# Transition

# Centroid R
hPiCentroidRTransition = TH1F('pi_centroid_r_transition', 'Pion Centroid Radius (Transition)', 500, 1800, 4200)
hists.append(hPiCentroidRTransition)

hElecCentroidRTransition = TH1F('elec_centroid_r_transition', 'Electron Centroid Radius (Transition)', 500, 1800, 4200)
hists.append(hElecCentroidRTransition)

# E/(E+H)
hPiEOverEPlusHTransition = TH1F('pi_e_over_e_plus_h_transition', 'Pion E/(E+H) (Transition)', 50, 0, 1.1)
hists.append(hPiEOverEPlusHTransition)

hElecEOverEPlusHTransition = TH1F('elec_e_over_e_plus_h_transition', 'Electron E/(E+H) (Transition)', 50, 0, 1.1)
hists.append(hElecEOverEPlusHTransition)

# Shower start layer
hPiShowerStartLayerTransition = TH1F('pi_shower_start_layer_transition', 'Pion Shower Start Layer (Transition)', 21, 0, 20)
hists.append(hPiShowerStartLayerTransition)

hElecShowerStartLayerTransition = TH1F('elec_shower_start_layer_transition', 'Electron Shower Start Layer (Transition)', 21, 0, 20)
hists.append(hElecShowerStartLayerTransition)

# Endcap

# Centroid R
hPiCentroidREndcap = TH1F('pi_centroid_r_endcap', 'Pion Centroid Radius (Endcap)', 500, 1800, 4200)
hists.append(hPiCentroidREndcap)

hElecCentroidREndcap = TH1F('elec_centroid_r_endcap', 'Electron Centroid Radius (Endcap)', 500, 1800, 4200)
hists.append(hElecCentroidREndcap)

# E/(E+H)
hPiEOverEPlusHEndcap = TH1F('pi_e_over_e_plus_h_endcap', 'Pion E/(E+H) (Endcap)', 50, 0, 1.1)
hists.append(hPiEOverEPlusHEndcap)

hElecEOverEPlusHEndcap = TH1F('elec_e_over_e_plus_h_endcap', 'Electron E/(E+H) (Endcap)', 50, 0, 1.1)
hists.append(hElecEOverEPlusHEndcap)

# Shower start layer
hPiShowerStartLayerEndcap = TH1F('pi_shower_start_layer_endcap', 'Pion Shower Start Layer (Endcap)', 21, 0, 20)
hists.append(hPiShowerStartLayerEndcap)

hElecShowerStartLayerEndcap = TH1F('elec_shower_start_layer_endcap', 'Electron Shower Start Layer (Endcap)', 21, 0, 20)
hists.append(hElecShowerStartLayerEndcap)

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
        clusters = event.getCollection('PandoraClusters')
        mc_particles = event.getCollection('MCParticle')

        for pfo in pfos:
            pfo_pdg = abs(pfo.getType())
            pfo_theta = getTheta(pfo)
            if (pfo_pdg == 11 or pfo_pdg == 211):
                matched_mc_particle = getMatchedMCParticle(pfo, mc_particles)
                if (matched_mc_particle):
                    matched_pdg = abs(matched_mc_particle.getPDG())
                    if ((matched_pdg == 211 or matched_pdg == 11) and fromTau(matched_mc_particle)):
                        clusters = pfo.getClusters()
                        if (len(clusters) > 0):
                            max_energy = 0
                            max_energy_cluster = None
                            for cluster in clusters:
                                energy = cluster.getEnergy()
                                if (energy > max_energy):
                                    max_energy = energy
                                    max_energy_cluster = cluster
                            centroid = max_energy_cluster.getPosition()
                            centroid_r = math.sqrt(centroid[0]**2 + centroid[1]**2 + centroid[2]**2)
                            ecal = max_energy_cluster.getSubdetectorEnergies()[0]
                            hcal = max_energy_cluster.getSubdetectorEnergies()[1]

                            if (matched_pdg == 211):
                                if (pfo_theta > 1 and pfo_theta < 2):
                                    hPiCentroidRCentralBarrel.Fill(centroid_r)
                                    hPiEOverEPlusHCentralBarrel.Fill(ecal/(ecal+hcal))
                                elif ((pfo_theta > 0.577 and pfo_theta < 1.0) or (pfo_theta > 2.0 and pfo_theta < 2.56)):
                                    hPiCentroidRTransition.Fill(centroid_r)
                                    hPiEOverEPlusHTransition.Fill(ecal/(ecal+hcal))
                                elif (pfo_theta < 0.577 or pfo_theta > 2.56):
                                    hPiCentroidREndcap.Fill(centroid_r)
                                    hPiEOverEPlusHEndcap.Fill(ecal/(ecal+hcal))
                            else:
                                if (pfo_theta > 1 and pfo_theta < 2):
                                    hElecCentroidRCentralBarrel.Fill(centroid_r)
                                    hElecEOverEPlusHCentralBarrel.Fill(ecal/(ecal+hcal))
                                elif ((pfo_theta > 0.577 and pfo_theta < 1.0) or (pfo_theta > 2.0 and pfo_theta < 2.56)):
                                    hElecCentroidRTransition.Fill(centroid_r)
                                    hElecEOverEPlusHTransition.Fill(ecal/(ecal+hcal))
                                elif (pfo_theta < 0.577 or pfo_theta > 2.56):
                                    hElecCentroidREndcap.Fill(centroid_r)
                                    hElecEOverEPlusHEndcap.Fill(ecal/(ecal+hcal))
                    
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
