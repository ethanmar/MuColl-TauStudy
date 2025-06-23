from pyLCIO import IOIMPL, EVENT, UTIL
from argparse import ArgumentParser
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

from tau_mc_link import getDecayMode, getLinkedMCTau

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('--tauInputFile', type=str, default='output_taufinder.slcio')
parser.add_argument('--pandoraInputFile', type=str, default='reco_output.slcio')

args = parser.parse_args()

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
