from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas
import math
from argparse import ArgumentParser
from array import array
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import matplotlib.patches as patches

from tau_mc_link import getLinkedMCTau, getDecayMode

def isTauDaughter(mc_pi):
    parents = mc_pi.getParents()
    for parent in parents:
        if (abs(parent.getPDG()) == 15):
            return True
    return False

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
parser.add_argument('--eventsFile', type=str, default='failed_events.txt')

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

with open(args.eventsFile, "r") as f:
    events = [int(line.strip()) for line in f if line.strip()]

# Open input file(s)
for file in to_process:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    # Loop through events
    for ievt, event in enumerate(reader):

        if event.getEventNumber() in events:
        
            # Get collections
            # clusters = event.getCollection('PandoraClusters')
            # tracks = event.getCollection('SiTracks_Refitted')
            # pfos = event.getCollection('PandoraPFOs')
            # mc_particles = event.getCollection('MCParticle')
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

            # Event dump
            '''print(f'\nEvent {event.getEventNumber()}')
            print(f'\nReconstructed particles:')
            for pfo in pfos:
                mom = pfo.getMomentum()
                pt = math.sqrt(mom[0]**2 + mom[1]**2)
                pfo_theta = math.atan2(pt, mom[2])
                pfo_phi = math.atan2(mom[1], mom[0])
                pfo_energy = pfo.getEnergy()
                pfo_type = pfo.getType()
                print(f'\nPFO: Type: {pfo_type}, Pt: {pt:.2f}, Energy: {pfo_energy:.2f}')
                clusters = pfo.getClusters()
                if len(clusters) > 0:
                    for cluster in clusters:
                        energy = cluster.getEnergy()
                        subdet_energies = cluster.getSubdetectorEnergies()
                        ecal_energy = subdet_energies[0]
                        hcal_energy = subdet_energies[1]
                        pos = cluster.getPosition()
                        r = math.sqrt(pos[0]**2 + pos[1]**2)
                        z = pos[2]
                        print(f'CLUSTER: Energy: {energy:.2f}, Ecal energy: {ecal_energy:.2f}, Hcal energy: {hcal_energy:.2f}, R: {r:.2f}, Z: {z:.2f}')
                tracks = pfo.getTracks()
                if len(tracks) > 0:
                    for track in tracks:
                        d0 = track.getD0()
                        omega = track.getOmega()
                        dEdx = track.getdEdx()
                        subdet_hits = track.getSubdetectorHitNumbers()
                        vtx_hits = subdet_hits[0]
                        in_hits = subdet_hits[1]
                        out_hits = subdet_hits[2]
                        print(f'TRACK: d0: {d0:.2f}, Omega: {omega:.5f}, dEdx: {dEdx:.5f}, Vertex hits: {vtx_hits}, Inner hits: {in_hits}, Outer hits: {out_hits}')
                for mc_particle in mc_particles:
                    pdg = mc_particle.getPDG()
                    mom = mc_particle.getMomentum()
                    pt = math.sqrt(mom[0]**2 + mom[1]**2)
                    mc_theta = math.atan2(pt, mom[2])
                    mc_phi = math.atan2(mom[1], mom[0])
                    deltaTheta = pfo_theta - mc_theta
                    deltaPhi = pfo_phi - mc_phi
                    while deltaPhi > math.pi:
                        deltaPhi -= 2*math.pi
                    while deltaPhi < -math.pi:
                        deltaPhi += 2*math.pi
                    deltaR = math.sqrt(deltaPhi**2 + deltaTheta**2)
                    mc_energy = mc_particle.getEnergy()
                    energy_diff = pfo_energy - mc_energy
                    endpoint = mc_particle.getEndpoint()
                    endpoint_r = math.sqrt(endpoint[0]**2 + endpoint[1]**2)
                    endpoint_z = endpoint[2]
                    if (deltaR < 0.1):
                        print(f'MCPARTICLE: PDG: {pdg}, Pt: {pt:.2f}, Energy: {mc_energy:.2f}, Endpoint r: {endpoint_r:.2f}, Endpoint z: {endpoint_z:.2f}, DeltaR: {deltaR:.5f}, Energy diff: {energy_diff:.2f}')

            print(f'\nMC Pions:')
            for mc_particle in mc_particles:
                if (abs(mc_particle.getPDG()) == 211 and isTauDaughter(mc_particle)):
                    mom = mc_particle.getMomentum()
                    pt = math.sqrt(mom[0]**2 + mom[1]**2)
                    energy = mc_particle.getEnergy()
                    endpoint = mc_particle.getEndpoint()
                    endpoint_r = math.sqrt(endpoint[0]**2 + endpoint[1]**2)
                    endpoint_z = endpoint[2]
                    print(f'\nMC PION: Pt: {pt:.2f}, Energy: {energy:.2f}, Endpoint r: {endpoint_r:.2f}, Endpoint z: {endpoint_z:.2f}')
                    daughters = mc_particle.getDaughters()
                    for daughter in daughters:
                        pdg = daughter.getPDG()
                        mom = daughter.getMomentum()
                        pt = math.sqrt(mom[0]**2 + mom[1]**2)
                        energy = daughter.getEnergy()
                        vertex = daughter.getVertex()
                        vertex_r = math.sqrt(vertex[0]**2 + vertex[1]**2)
                        vertex_z = vertex[2]
                        vertex_not_endpoint = (1 if daughter.vertexIsNotEndpointOfParent() else 0)
                        print(f'DAUGHTER: PDG: {pdg}, Pt: {pt:.2f}, Energy: {energy:.2f}, Vertex r: {vertex_r:.2f}, Vertex z: {vertex_z:.2f}, Vertex not endpoint: {vertex_not_endpoint}')'''
                        
            # Event display
            E_y_x_digi = np.zeros((len(y_bins)-1, len(x_bins)-1))
            E_r_z_digi = np.zeros((len(r_bins)-1, len(z_bins)-1))
                
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
                E_y_x_digi[y_bin][x_bin] += E
                E_r_z_digi[r_bin][z_bin] += E
                
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
                E_y_x_digi[y_bin][x_bin] += E
                E_r_z_digi[r_bin][z_bin] += E

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
                E_y_x_digi[y_bin][x_bin] += E
                E_r_z_digi[r_bin][z_bin] += E

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
                E_y_x_digi[y_bin][x_bin] += E
                E_r_z_digi[r_bin][z_bin] += E

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
                E_y_x_digi[y_bin][x_bin] += E
                E_r_z_digi[r_bin][z_bin] += E

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
                E_y_x_digi[y_bin][x_bin] += E
                E_r_z_digi[r_bin][z_bin] += E

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
                    E_y_x_digi[y_bin][x_bin] += E
                    E_r_z_digi[r_bin][z_bin] += E

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
                    E_y_x_digi[y_bin][x_bin] += E
                    E_r_z_digi[r_bin][z_bin] += E

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
                    E_y_x_digi[y_bin][x_bin] += E
                    E_r_z_digi[r_bin][z_bin] += E

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
                    E_y_x_digi[y_bin][x_bin] += E
                    E_r_z_digi[r_bin][z_bin] += E

            masked_E_y_x_digi = np.ma.masked_where(E_y_x_digi == 0, E_y_x_digi)
            masked_E_r_z_digi = np.ma.masked_where(E_r_z_digi == 0, E_r_z_digi)            

            plt.close()

            plt.figure()
            cmap = plt.cm.jet.copy()
            cmap.set_bad(color='white')

            plt.clf()
            plt.imshow(masked_E_y_x_digi, extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]], aspect='auto', origin='lower', cmap=cmap)
            plt.xticks(np.arange(x_bins[0], x_bins[-1]+1, 1000), rotation=90)
            plt.yticks(np.arange(y_bins[0], y_bins[-1]+1, 1000))
            ax = plt.gca()
            ax.xaxis.set_minor_locator(MultipleLocator(200))
            ax.yaxis.set_minor_locator(MultipleLocator(200))
            ax.tick_params(axis='both', which='major', length=7, labelsize=10)
            ax.tick_params(axis='both', which='minor', length=3)

            ax = plt.gca()
            solenoid = patches.Wedge((0, 0), 1857, 0, 360, width=1857-1500, edgecolor='violet', facecolor='violet', alpha=0.2, label='Solenoid')
            ax.add_patch(solenoid)
            ecal = patches.Wedge((0, 0), 2125, 0, 360, width=2125-1857, edgecolor='deeppink', facecolor='deeppink', alpha=0.2, label='ECal')
            ax.add_patch(ecal)
            hcal = patches.Wedge((0, 0), 4113, 0, 360, width=4113-2125, edgecolor='darkviolet', facecolor='darkviolet', alpha=0.2, label='HCal')
            ax.add_patch(hcal)

            plt.xlabel('X [mm]')
            plt.ylabel('Y [mm]')
            plt.colorbar(label='Energy [GeV]')
            plt.legend(loc='upper right')
            plt.tight_layout()
            plt.savefig(f'x_y_E_digi_{event.getEventNumber()}.png')

            plt.clf()
            plt.imshow(masked_E_r_z_digi, extent=[z_bins[0], z_bins[-1], r_bins[0], r_bins[-1]], aspect='auto', origin='lower', cmap=cmap)
            plt.xticks(np.arange(z_bins[0], z_bins[-1]+1, 1000), rotation=90)
            plt.yticks(np.arange(r_bins[0], r_bins[-1]+1, 500))
            ax = plt.gca()
            ax.xaxis.set_minor_locator(MultipleLocator(200))
            ax.yaxis.set_minor_locator(MultipleLocator(100))
            ax.tick_params(axis='both', which='major', length=7, labelsize=10)
            ax.tick_params(axis='both', which='minor', length=3)

            ax = plt.gca()
            solenoid = patches.Rectangle((-2307, 1500), 2*2307, 1857-1500, edgecolor='violet', facecolor='violet', alpha=0.2, label='Solenoid')
            ax.add_patch(solenoid)
            ecal_barrel = patches.Rectangle((-2307, 1857), 2*2307, 2125-1857, edgecolor='deeppink', facecolor='deeppink', alpha=0.2, label='ECal')
            ax.add_patch(ecal_barrel)
            ecal_endcap_left = patches.Rectangle((-2575, 310), 2575-2307, 2125-310, edgecolor='deeppink', facecolor='deeppink', alpha=0.2)
            ax.add_patch(ecal_endcap_left)
            ecal_endcap_right = patches.Rectangle((2307, 310), 2575-2307, 2125-310, edgecolor='deeppink', facecolor='deeppink', alpha=0.2)
            ax.add_patch(ecal_endcap_right)
            hcal_barrel = patches.Rectangle((-2575, 2125), 2*2575, 4113-2125, edgecolor='darkviolet', facecolor='darkviolet', alpha=0.2, label='HCal')
            ax.add_patch(hcal_barrel)
            hcal_endcap_left = patches.Rectangle((-4562, 307), 4562-2575, 4112-307, edgecolor='darkviolet', facecolor='darkviolet', alpha=0.2)
            ax.add_patch(hcal_endcap_left)
            hcal_endcap_right = patches.Rectangle((2575, 307), 4562-2575, 4112-307, edgecolor='darkviolet', facecolor='darkviolet', alpha=0.2)
            ax.add_patch(hcal_endcap_right)
            
            plt.xlabel('Z [mm]')
            plt.ylabel('R [mm]')
            plt.colorbar(label='Energy [GeV]')
            plt.legend(loc='upper right')
            plt.tight_layout()
            plt.savefig(f'z_r_E_digi_{event.getEventNumber()}.png')
            
    reader.close()
