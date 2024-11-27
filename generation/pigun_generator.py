#####################################
#
# simple script to create lcio files with single particle
# events - modify as needed
# @author F.Gaede, DESY
# @date 1/07/2014
#
# initialize environment:
#  export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib
#
#####################################
import numpy as np
import random
from array import array
from g4units import deg, s

# --- LCIO dependencies ---
from pyLCIO import EVENT, IMPL, IOIMPL

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--numberOfEvents", type=int, default=1000)
# parser.add_argument("--gunEnergy", type = float)
parser.add_argument("--output", type=str, default="piguns.slcio")
args = parser.parse_args()

# ---- number of events per momentum bin -----
nevt = args.numberOfEvents

outfile = args.output

# --------------------------------------------

wrt = IOIMPL.LCFactory.getInstance().createLCWriter()

wrt.open(outfile, EVENT.LCIO.WRITE_NEW)

random.seed()


# ========== particle properties ===================

energies = []

genstat = 1 #change based on your configuration
pdg = 211
mass = 0.13957  # pion mass
charge = +1

# decay time in seconds
lifetime = 2.6033e-8 * s

# bounds on theta
theta_min = 10.0 * deg
theta_max = 170.0 * deg

# =================================================

for j in range(0, nevt):
    col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE)
    evt = IMPL.LCEventImpl()

    evt.setEventNumber(j)

    evt.addCollection(col, "MCParticle")
    
    # --------- generate particle properties ----------
    
    E = random.uniform(5., 300.) #flat in E between 5 and 300 GeV

    phi = random.random() * np.pi * 2.0  # flat in phi
    
    theta = random.uniform(theta_min, theta_max) # flat in theta

    p = np.sqrt(E**2 - mass**2)

    px = p * np.sin(theta) * np.cos(phi)
    py = p * np.sin(theta) * np.sin(phi)
    pz = p * np.cos(theta)

    momentum = array("f", [px, py, pz])

    # --------------- create MCParticle -------------------

    mcp = IMPL.MCParticleImpl()

    mcp.setGeneratorStatus(genstat)
    mcp.setMass(mass)
    mcp.setPDG(pdg)
    mcp.setMomentum(momentum)
    mcp.setCharge(charge)

    # -------------------------------------------------------

    col.addElement(mcp)

    wrt.writeEvent(evt)


wrt.close()
