'''
This script reads Cherenkov and fluorescence photon bunches
from all events in the CORSIKA output.

Usage: It takes CERnnnnnn and DATACARD files as input arguments.

    python test.py CERnnnnnn all

Variables in capital letters corespond to datacard input arguments
'''

from corsika import histogram, DataCard
import sys
import os
import argparse
from itertools import compress
from math import pi
import numpy as np
from scipy.io import FortranFile


data_card = DataCard.read(sys.argv[2])
# Telescope pointing angle
onaxis = input("On-axis pointing (y/n)? ")
if onaxis == 'y':
    pointing_angle = data_card['THETAP']
    pointing = 'onaxis'
else:
    pointing_angle = input("Off-axis angle? ")
    pointing = 'offaxis'

# Histogram definition depending on the detection area
binsize = 10 #meters

if data_card['XCERARY'] > data_card['XCERARY']:
    print("Histogram along x-axis...")
    type_of_hist = 'x'
    maxlen = 1e-2 * data_card['XCERARY'] / 2
    numbins = int(maxlen / binsize)
    breaks = numbins + 1
    distances = np.linspace(0, maxlen, breaks)
    mids = ((distances[1] - distances[0]) / 2) + distances
    mids = mids[0:numbins]

elif data_card['XCERARY'] < data_card['YCERARY']:
    print("Histogram along y-axis...")
    type_of_hist = 'y'
    maxlen = 1e-2 * data_card['YCERARY'] / 2
    numbins = int(maxlen / binsize)
    breaks = numbins + 1
    distances = np.linspace(0, maxlen, breaks)
    mids = ((distances[1] - distances[0]) / 2) + distances
    mids = mids[0:numbins]

else:
    print("Histogram radially...")
    type_of_hist = 'r'
    maxlen  = 1e-2 * data_card['XCERARY'] / 2
    numbins = int(maxlen / binsize)
    breaks = numbins + 1
    radius = np.linspace(0, maxlen, breaks)
    mids = ((radius[1] - radius[0]) / 2) + radius
    mids = mids[0:numbins]
    ring = pi * ((radius[1:]) ** 2 - (radius[0:numbins]) ** 2)

bunches = np.array([]).reshape(0, 7)
hist_c = np.zeros((2, numbins))
hist_f = np.zeros((2, numbins))

# Open binary CORSIKA output file
file = FortranFile(sys.argv[1], 'r')

# Control and debugging counters:
count = 0
lines = 0
photons = 0

while True:
    count = count + 1
    # Sort data in 21 sub-blocks of 39 lines each (and 7 columns)
    data = np.split(file.read_reals(dtype=np.float32).reshape(-1, 7), 21)
    # It should be:
    # indices_boolean = [np.abs(i[0][0])<max(cersize, fluorsize) for i in data]
    # select only sub-blocks of bunches
    indices_boolean = [np.abs(i[0][0]) < 100 for i in data]
    # store 1st element of each sub-block
    indices = [i[0][0] for i in data]
    bunches = np.vstack([bunches, np.vstack(compress(data, indices_boolean))])
    # drop those lines containing only zeros
    bunches = bunches[np.all(bunches != 0, axis=1)]

    # Not store sub-blocks to bunches array every time to speed the process up
    if count == 10:
        h_c = histogram(bunches[bunches[:, 0] > 0]
                        , data_card['XCERARY']
                        , data_card['YCERARY']
                        , pointing_angle)
        h_f = histogram(bunches[bunches[:, 0] < 0]
                        , data_card['XCERARY']
                        , data_card['YCERARY']
                        , pointing_angle)
        hist_c = hist_c + h_c
        hist_f = hist_f + h_f
        count = 0  # reset counter
        bunches = np.array([]).reshape(0, 7)  # reset array

    if any(3300 < i < 3303. for i in indices):  # Flag indicating RUN END sub-block
        if count < 10:
            h_c = histogram(bunches[bunches[:, 0] > 0]
                            , data_card['XCERARY']
                            , data_card['YCERARY']
                            , pointing_angle)
            h_f = histogram(bunches[bunches[:, 0] < 0]
                            , data_card['XCERARY']
                            , data_card['YCERARY']
                            , pointing_angle)
            hist_c = hist_c + h_c
            hist_f = hist_f + h_f
        file.close()
        break

# if __name__ == '__main__':
#     # parse command line arguments
#     parser = argparse.ArgumentParser(description='Read photon bunches from CORSIKA standard output file')
#     parser.add_argument('CERnnnnnn', help='path to CERnnnnnn binary file containing data')
#     parser.add_argument('--datacard', default='all')
#     args = parser.parse_args()
#
#     # train(variable_input_model, args.h5_file, args.epochs, args.image_summary, args.embedding)
