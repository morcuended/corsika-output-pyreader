'''
This script reads the input arguments passed in the DATA CARD.

'''

import numpy as np
from astropy.io import ascii
import os, sys

# =========  DATA CARD =========
with open(sys.argv[2], 'r') as f:
    read_data = f.read()
    #print(read_data) 
    datacard = ascii.read(read_data, format='fixed_width_no_header'
                          , delimiter=' '
                          , col_starts=(0, 8, 39))
f.closed

nshower = int(datacard[2][1].split()[0])
E_prim = float(datacard[5][1].split()[1])
prim_part = int(datacard[3][1].split()[0])
if prim_part == 1:
    prim_part = 'gamma'
if prim_part == 14:
    prim_part = 'proton'
seed1 = int(datacard[8][1].split()[0])
seed2 = int(datacard[9][1].split()[0])
theta_p = float(datacard[6][1].split()[0])
phi_p = float(datacard[7][1].split()[0])
obs_level = float(datacard[10][1].split()[0])
obs_level = obs_level * 1e-2 #in meters
atm_mod = int(datacard[18][1].split()[0])
cersize = float(datacard[23][1].split()[0])
fluorsize = float(datacard[24][1].split()[0])
x_area = float(datacard[27][1].split()[4])
y_area = float(datacard[27][1].split()[5])


