from tau_mc_link import getDecayMode

from pyLCIO import IOIMPL, EVENT, UTIL
from argparse import ArgumentParser
import os

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--inputFile', type=str, default='output_taufinder.slcio')

args = parser.parse_args()

# Create list to store frequency of decay modes
decay_modes = [0, 0, 0, 0, 0, 0, 0, 0]

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

        #Get MCParticle collection
        mcParticles = event.getCollection('MCParticle')

        # Loop through mcParticles
        for mcParticle in mcParticles:

            # Get mcTaus
            if abs(mcParticle.getPDG()) == 15:

                # Get decay mode
                decay_mode = getDecayMode(mcParticle)
                if decay_mode is not None:
                    decay_modes[decay_mode] += 1

    # Close file
    reader.close()

# Open and write txt file
with open('branching.txt', 'w') as f:
    f.write(f'pi-, nu-tau: {(decay_modes[0]/sum(decay_modes))*100}%\n')
    f.write(f'pi-, pi0, nu-tau: {(decay_modes[1]/sum(decay_modes))*100}%\n')
    f.write(f'pi-, pi0x2, nu-tau: {(decay_modes[2]/sum(decay_modes))*100}%\n')
    f.write(f'pi-, pi0x3, nu-tau: {(decay_modes[3]/sum(decay_modes))*100}%\n')
    f.write(f'3-prong plus nu-tau: {(decay_modes[4]/sum(decay_modes))*100}%\n')
    f.write(f'3-prong plus pi0 and nu-tau: {(decay_modes[5]/sum(decay_modes))*100}%\n')
    f.write(f'nu-e-bar, e-, nu-tau: {(decay_modes[6]/sum(decay_modes))*100}%\n')
    f.write(f'nu-mu-bar, mu-, nu-tau: {(decay_modes[7]/sum(decay_modes))*100}%\n')


