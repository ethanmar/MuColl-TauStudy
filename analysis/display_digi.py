from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

from tau_mc_link import getLinkedMCTau, getDecayMode

def getXBin(x, x_bins):
    for i in range(len(x_bins)-1):
        if (x >= x_bins[i] and x < x_bins[i+1]):
            return i
    return len(x_bins)-2

def getYBin(y, y_bins):
    for i in range(len(y_bins)-1):
        if (y >= y_bins[i] and y < y_bins[i+1]):
            return i
    return len(y_bins)-2

def getZBin(z, z_bins):
    for i in range(len(z_bins)-1):
        if (z >= z_bins[i] and z < z_bins[i+1]):
            return i
    return len(z_bins)-2

def getRBin(r, r_bins):
    for i in range(len(r_bins)-1):
        if (r >= r_bins[i] and r < r_bins[i+1]):
            return i
    return len(r_bins)-2

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_reco.slcio')

args = parser.parse_args()

# Check if input file is a directory or a single file    
to_process = []

if os.path.isdir(args.inputFile):
    for r, d, f in os.walk(args.inputFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(args.inputFile)

x_bins = np.linspace(-5000, 5000, 101)
y_bins = np.linspace(-5000, 5000, 101)
z_bins = np.linspace(-5000, 5000, 51)
r_bins = np.linspace(0, 5000, 51)


events = [354, 373, 468, 475, 663, 685, 730, 1796, 2219, 2454, 2885, 2944, 3023, 3111, 3275, 3497, 3706, 3982, 11434, 11935, 12178, 12527, 13128, 13218, 13489, 13540, 13924, 13988, 14229, 14325, 14463, 14605, 14657, 3999, 4211, 4526, 4568, 4781, 5383, 5551, 6340, 7035, 7142, 7166, 7374, 7459, 7764, 7825, 7918, 7965, 7966, 7994, 8038, 8165, 8428, 8529, 8639, 8651, 8760, 9431, 9470, 9558, 9837, 10090, 10399, 10651, 10664, 10940, 11389]

# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        if event.getEventNumber() in events:
        
            # Get collections
            vertex_barrel_hits = event.getCollection('VBTrackerHitsConed')
            vertex_endcap_hits = event.getCollection('VETrackerHitsConed')
            inner_trkr_barrel_hits = event.getCollection('IBTrackerHitsConed')
            inner_trkr_endcap_hits = event.getCollection('IETrackerHitsConed')
            outer_trkr_barrel_hits = event.getCollection('OBTrackerHitsConed')
            outer_trkr_endcap_hits = event.getCollection('OETrackerHitsConed')
            if 'EcalBarrelCollectionSel' in event.getCollectionNames():
                ecal_barrel_hits = event.getCollection('EcalBarrelCollectionSel')
            else:
                ecal_barrel_hits = None
            if 'EcalEndcapCollectionSel' in event.getCollectionNames():
                ecal_endcap_hits = event.getCollection('EcalEndcapCollectionSel')
            else:
                ecal_endcap_hits = None
            if 'HcalBarrelCollectionConed' in event.getCollectionNames():
                hcal_barrel_hits = event.getCollection('HcalBarrelCollectionConed')
            else:
                hcal_barrel_hits = None
            if 'HcalEndcapCollectionConed' in event.getCollectionNames():
                hcal_endcap_hits = event.getCollection('HcalEndcapCollectionConed')
            else:
                hcal_endcap_hits = None

            E_y_x = np.zeros((len(y_bins)-1, len(x_bins)-1))
            E_r_z = np.zeros((len(r_bins)-1, len(z_bins)-1))
                
            # Loop through vertex barrel hits
            for vertex_barrel_hit in vertex_barrel_hits:
                # Get hit position
                pos = vertex_barrel_hit.getPosition()
                x = pos[0]
                x_bin = getXBin(x, x_bins)
                y = pos[1]
                y_bin = getYBin(y, y_bins)
                z = pos[2]
                z_bin = getZBin(z, z_bins)
                r = math.sqrt(x**2 + y**2 + z**2)
                r_bin = getRBin(r, r_bins)
                E = vertex_barrel_hit.getEDep()
                E_y_x[y_bin][x_bin] += E
                E_r_z[r_bin][z_bin] += E
                
            # Loop through vertex endcap hits
            for vertex_endcap_hit in vertex_endcap_hits:
                # Get hit position
                pos = vertex_endcap_hit.getPosition()
                x = pos[0]
                x_bin = getXBin(x, x_bins)
                y = pos[1]
                y_bin = getYBin(y, y_bins)
                z = pos[2]
                z_bin = getZBin(z, z_bins)
                r = math.sqrt(x**2 + y**2 + z**2)
                r_bin = getRBin(r, r_bins)
                E = vertex_endcap_hit.getEDep()
                E_y_x[y_bin][x_bin] += E
                E_r_z[r_bin][z_bin] += E

            # Loop through inner tracker barrel hits
            for inner_trkr_barrel_hit in inner_trkr_barrel_hits:
                # Get hit position
                pos = inner_trkr_barrel_hit.getPosition()
                x = pos[0]
                x_bin = getXBin(x, x_bins)
                y = pos[1]
                y_bin = getYBin(y, y_bins)
                z = pos[2]
                z_bin = getZBin(z, z_bins)
                r = math.sqrt(x**2 + y**2 + z**2)
                r_bin = getRBin(r, r_bins)
                E = inner_trkr_barrel_hit.getEDep()
                E_y_x[y_bin][x_bin] += E
                E_r_z[r_bin][z_bin] += E

            # Loop through inner tracker endcap hits
            for inner_trkr_endcap_hit in inner_trkr_endcap_hits:
                # Get hit position
                pos = inner_trkr_endcap_hit.getPosition()
                x = pos[0]
                x_bin = getXBin(x, x_bins)
                y = pos[1]
                y_bin = getYBin(y, y_bins)
                z = pos[2]
                z_bin = getZBin(z, z_bins)
                r = math.sqrt(x**2 + y**2 + z**2)
                r_bin = getRBin(r, r_bins)
                E = inner_trkr_endcap_hit.getEDep()
                E_y_x[y_bin][x_bin] += E
                E_r_z[r_bin][z_bin] += E

            # Loop through outer tracker barrel hits
            for outer_trkr_barrel_hit in outer_trkr_barrel_hits:
                # Get hit position
                pos = outer_trkr_barrel_hit.getPosition()
                x = pos[0]
                x_bin = getXBin(x, x_bins)
                y = pos[1]
                y_bin = getYBin(y, y_bins)
                z = pos[2]
                z_bin = getZBin(z, z_bins)
                r = math.sqrt(x**2 + y**2 + z**2)
                r_bin = getRBin(r, r_bins)
                E = outer_trkr_barrel_hit.getEDep()
                E_y_x[y_bin][x_bin] += E
                E_r_z[r_bin][z_bin] += E

            # Loop through outer tracker endcap hits
            for outer_trkr_endcap_hit in outer_trkr_endcap_hits:
                # Get hit position
                pos = outer_trkr_endcap_hit.getPosition()
                x = pos[0]
                x_bin = getXBin(x, x_bins)
                y = pos[1]
                y_bin = getYBin(y, y_bins)
                z = pos[2]
                z_bin = getZBin(z, z_bins)
                r = math.sqrt(x**2 + y**2 + z**2)
                r_bin = getRBin(r, r_bins)
                E = outer_trkr_endcap_hit.getEDep()
                E_y_x[y_bin][x_bin] += E
                E_r_z[r_bin][z_bin] += E

            # Loop through ecal barrel hits
            if ecal_barrel_hits is not None:
                for ecal_barrel_hit in ecal_barrel_hits:
                    # Get hit position
                    pos = ecal_barrel_hit.getPosition()
                    x = pos[0]
                    x_bin = getXBin(x, x_bins)
                    y = pos[1]
                    y_bin = getYBin(y, y_bins)
                    z = pos[2]
                    z_bin = getZBin(z, z_bins)
                    r = math.sqrt(x**2 + y**2 + z**2)
                    r_bin = getRBin(r, r_bins)
                    E = ecal_barrel_hit.getEnergy()
                    E_y_x[y_bin][x_bin] += E
                    E_r_z[r_bin][z_bin] += E

            # Loop through ecal endcap hits
            if ecal_endcap_hits is not None:
                for ecal_endcap_hit in ecal_endcap_hits:
                    # Get hit position
                    pos = ecal_endcap_hit.getPosition()
                    x = pos[0]
                    x_bin = getXBin(x, x_bins)
                    y = pos[1]
                    y_bin = getYBin(y, y_bins)
                    z = pos[2]
                    z_bin = getZBin(z, z_bins)
                    r = math.sqrt(x**2 + y**2 + z**2)
                    r_bin = getRBin(r, r_bins)
                    E = ecal_endcap_hit.getEnergy()
                    E_y_x[y_bin][x_bin] += E
                    E_r_z[r_bin][z_bin] += E

            # Loop through hcal barrel hits
            if hcal_barrel_hits is not None:
                for hcal_barrel_hit in hcal_barrel_hits:
                    # Get hit position
                    pos = hcal_barrel_hit.getPosition()
                    x = pos[0]
                    x_bin = getXBin(x, x_bins)
                    y = pos[1]
                    y_bin = getYBin(y, y_bins)
                    z = pos[2]
                    z_bin = getZBin(z, z_bins)
                    r = math.sqrt(x**2 + y**2 + z**2)
                    r_bin = getRBin(r, r_bins)
                    E = hcal_barrel_hit.getEnergy()
                    E_y_x[y_bin][x_bin] += E
                    E_r_z[r_bin][z_bin] += E

            # Loop through hcal endcap hits
            if hcal_endcap_hits is not None:
                for hcal_endcap_hit in hcal_endcap_hits:
                    # Get hit position
                    pos = hcal_endcap_hit.getPosition()
                    x = pos[0]
                    x_bin = getXBin(x, x_bins)
                    y = pos[1]
                    y_bin = getYBin(y, y_bins)
                    z = pos[2]
                    z_bin = getZBin(z, z_bins)
                    r = math.sqrt(x**2 + y**2 + z**2)
                    r_bin = getRBin(r, r_bins)
                    E = hcal_endcap_hit.getEnergy()
                    E_y_x[y_bin][x_bin] += E
                    E_r_z[r_bin][z_bin] += E

            masked_E_y_x = np.ma.masked_where(E_y_x == 0, E_y_x)
            masked_E_r_z = np.ma.masked_where(E_r_z == 0, E_r_z)            

            plt.figure()
            cmap = plt.cm.jet.copy()
            cmap.set_bad(color='white')

            plt.clf()
            plt.imshow(masked_E_y_x, extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]], aspect='auto', origin='lower', cmap=cmap)
            plt.xticks(np.arange(x_bins[0], x_bins[-1]+1, 1000))
            plt.yticks(np.arange(y_bins[0], y_bins[-1]+1, 1000))
            ax = plt.gca()
            ax.xaxis.set_minor_locator(MultipleLocator(200))
            ax.yaxis.set_minor_locator(MultipleLocator(200))
            ax.tick_params(axis='both', which='major', length=7, labelsize=10)
            ax.tick_params(axis='both', which='minor', length=3)
            plt.xlabel('X [mm]')
            plt.ylabel('Y [mm]')
            plt.colorbar(label='E [GeV]')
            plt.tight_layout()
            plt.savefig(f'x_y_E_{event.getEventNumber()}.png')

            plt.clf()
            plt.imshow(masked_E_r_z, extent=[z_bins[0], z_bins[-1], r_bins[0], r_bins[-1]], aspect='auto', origin='lower', cmap=cmap)
            plt.xticks(np.arange(z_bins[0], z_bins[-1]+1, 1000))
            plt.yticks(np.arange(r_bins[0], r_bins[-1]+1, 500))
            ax = plt.gca()
            ax.xaxis.set_minor_locator(MultipleLocator(200))
            ax.yaxis.set_minor_locator(MultipleLocator(100))
            ax.tick_params(axis='both', which='major', length=7, labelsize=10)
            ax.tick_params(axis='both', which='minor', length=3)
            plt.xlabel('Z [mm]')
            plt.ylabel('R [mm]')
            plt.colorbar(label='E [GeV]')
            plt.tight_layout()
            plt.savefig(f'z_r_E_{event.getEventNumber()}.png')

    reader.close()
