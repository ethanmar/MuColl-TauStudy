# Create lcio files with single tau events
# E. Martinez, 10/01/2024

# initialize environment:
# export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib

import math
import random
from array import array
from g4units import deg

# LCIO dependencies
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL

# Number of events per momentum bin
nevt = 1000

# Output file
outfile = "tau_gen.slcio"

# Open output file
wrt = IOIMPL.LCFactory.getInstance().createLCWriter( )
wrt.open(outfile, EVENT.LCIO.WRITE_NEW) 

#========== particle properties ===================
# random.seed()

'''
#generate uniform distribution in pT
pT = []
for i in range(500):
    pT.append(random.random()*300.+20) # 20-320 GeV/c
'''
pT = 100 # GeV/c
genstat  = 1
pdg = 15 # tau
mass =  1.77686 # GeV/c^2
charge = -1.
decayLen = 1.e22 
#=================================================



#for pt in pT:
for n in range(nevt):
    col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE) 
    evt = IMPL.LCEventImpl() 

    evt.setEventNumber(n)
    evt.addCollection(col, "MCParticle")

    phi =  random.random()*math.pi*2. #uniform distribution in phi
    theta = math.acos(2*random.random()-1) #uniform distribution in theta
                
    px = pT*math.cos(phi)
    py = pT*math.sin(phi)
    pz = pT/math.tan(theta) 

    # Momentum vector
    momentum  = array('f', [ px, py, pz ])  

    epx = decayLen*math.cos(phi)*math.sin(theta) 
    epy = decayLen*math.sin(phi)*math.sin(theta)
    epz = decayLen*math.cos(theta) 

    # Decay position vector
    endpoint = array('d',[ epx, epy, epz ] )  
        

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

    col.addElement( mcp )
        
    wrt.writeEvent( evt )


wrt.close() 
