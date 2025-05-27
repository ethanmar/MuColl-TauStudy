# Create lcio file with single tau events
# E. Martinez, 10/01/2024

import math
import numpy as np
import random
from array import array
from g4units import deg
import argparse

# LCIO dependencies
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL

# Command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--nEvents", type=int, default=1000)
parser.add_argument("--outputFile", type=str, default="tau_gen.slcio")
args = parser.parse_args()

# Number of events per momentum bin
nevt = args.nEvents

# Output file
outfile = args.outputFile

# Open output file
wrt = IOIMPL.LCFactory.getInstance().createLCWriter( )
wrt.open(outfile, EVENT.LCIO.WRITE_NEW) 

#========== particle properties ===================

# Set random seed
random.seed()

# Generator status
genstat  = 1

# PDG
pdg = 15 # tau

# Mass
mass =  1.77686 # GeV/c^2

# Charge
charge = -1.

# Decay length
decayLen = 1.e22

# Bounds on theta
theta_min = 10.0 * deg
theta_max = 170.0 * deg

#=================================================

# Loop over events
for n in range(nevt):
    
    # Initialize MCParticle collection and event
    col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE) 
    evt = IMPL.LCEventImpl() 

    # Set event number and add MCParticle collection
    evt.setEventNumber(n)
    evt.addCollection(col, "MCParticle")

    # Generate particle properties
    pT = random.random()*300.+20 # uniform in pT (20-320 GeV/c)

    phi =  random.random()*np.pi*2. # uniform in phi
    theta = random.uniform(theta_min, theta_max) # uniform in theta
                
    px = pT*np.cos(phi)
    py = pT*np.sin(phi)
    pz = pT/np.tan(theta) 

    momentum  = array('f', [px, py, pz])  

    epx = decayLen*np.cos(phi)*np.sin(theta) 
    epy = decayLen*np.sin(phi)*np.sin(theta)
    epz = decayLen*np.cos(theta) 

    # Decay position vector
    endpoint = array('d', [epx, epy, epz])  
        

#--------------- create MCParticle -------------------

    mcp = IMPL.MCParticleImpl() 

    mcp.setGeneratorStatus(genstat) 
    mcp.setMass(mass)
    mcp.setPDG(pdg) 
    mcp.setMomentum(momentum)
    mcp.setCharge(charge) 

    if(decayLen < 1.e9):   
        mcp.setEndpoint(endpoint)
        
#-------------------------------------------------------

    col.addElement(mcp)
        
    wrt.writeEvent(evt)


wrt.close() 
