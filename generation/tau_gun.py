# Create lcio file with single tau events
# E. Martinez, 10/01/2024

import math
import random
from array import array

# LCIO dependencies
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL

# Output file
outfile = "tau_gen_50k.slcio"

# Open output file
wrt = IOIMPL.LCFactory.getInstance().createLCWriter( )
wrt.open(outfile, EVENT.LCIO.WRITE_NEW) 

# Number of events
nevt = 50000

#========== particle properties ===================

random.seed()

genstat  = 1

pt_min = 20
pt_max = 320

theta_min = 8./180.*math.pi
theta_max = 172./180.*math.pi

pdg = 15 

mass =  1.77686 

charge = -1.

decayLen = 1.e22

#=================================================

for n in range(nevt):

    col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE)
    evt = IMPL.LCEventImpl()

    evt.setEventNumber(n)
    evt.addCollection(col, "MCParticle")


    pt = random.uniform(pt_min, pt_max)
    theta = random.uniform(theta_min, theta_max)
    phi = random.random() * math.pi * 2.
    p = pt/math.sin(theta)
    px = pt * math.cos(phi)
    py = pt * math.sin(phi)
    pz = p * math.cos(theta)

    momentum = array('f', [px, py, pz])

    epx = decayLen*math.cos(phi)*math.sin(theta)
    epy = decayLen*math.sin(phi)*math.sin(theta)
    epz = decayLen*math.cos(theta)

    endpoint = array('d', [epx, epy, epz])

# --------------- create MCParticle -------------------

    mcp = IMPL.MCParticleImpl()

    mcp.setGeneratorStatus(genstat)
    mcp.setMass(mass)
    mcp.setPDG(pdg)
    mcp.setMomentum(momentum)
    mcp.setCharge(charge)
        
    if (decayLen < 1.e9):   # arbitrary ...
        mcp.setEndpoint(endpoint)

# -------------------------------------------------------

    col.addElement(mcp)

    wrt.writeEvent(evt)

wrt.close()
