from scipy.io import FortranFile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi
from itertools import compress
from astropy.io import ascii
import matplotlib.style
import matplotlib as mpl
import sys

''' 
Ways to histogram showers (analogous to indicate where the telescope
is pointing):
  - on axis if theta_telesc = theta_primary
  - off axis (else)
  
Run this script as:

    python reader.py <filename_from_CORSIKA> <DATACARDname>

Then sys.argv[1] gets <filename_from_CORSIKA> as the 
file that will be analyze here and so sys.argv[2] gets the DATACARD
Also, the "all" DATA CARD file must be in the same directory.

TO DO list:
  (x) Establish a way of pointing the telescope to any direction
      via theta angle. 
	(x) Ask if user wants onaxis pointing, if not, ask for the pointing 
	  angle of the telescope.
  (-) Convert pointing angle to off-axis angle with respect to the shower axis direction.
  (-) If any argument is passed, CER000001 should be set
      as default input file.
  (x) Distinguish fluor/Cherenkov photons in order to histogramming both 
      components separately.
  (-) Store photon bunches w/ and w/o FoV constraint at the same time? 
      Would this be useful/efficient in the future (if we only want to store FoV-constrained photons)?
  (-) Average x/y-hitogramming.The idea is to have to diferent DATACARD with detector area defined 
      along both axis (x/y) with the corresponding shower direction for each one.
  (-) Check which dimension is the largest one, then histogramming along
      that direction.
'''
print('')
print('-------  CORSIKA reader  -------')

#with_fov = input("Include FoV (y/n)? ")

# Function to histogram Cherenkov photons
def histogram(bunches, x_area, y_area, theta):
    fov = np.cos(5 * pi / 180) #FoV contraint (+/- 5 deg)
    # Histogramming along x-axis
    if x_area > y_area:  
        weighted_pht = np.abs(bunches[:, 0]) / (binsize**2 * nshower)
#        if with_fov == 'n': # all photons
        h,edges = np.histogram(1e-2 * bunches[:, 1]
	              	       , bins = distances
	   		       , weights = weighted_pht
	   		       , range = [0., maxlen]
                               )
#        else: # 10 deg FoV
        wemis = np.sqrt(1-bunches[:, 3]**2 - bunches[:, 4]**2)
        # Pointing the telescope along the shower direction (on-axis) 
        # if theta(tel)!=theta_p -> off-axis observation
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180) 
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        weighted_pht_fov = weighted_pht[wemis >= fov]
        h_fov,edges = np.histogram(1e-2 * bunches_fov[:, 1]
	    		           , bins = distances
	    		           , weights = weighted_pht_fov
	    		           , range = [0., maxlen]
                                   )
    elif x_area < y_area:  
        weighted_pht = np.abs(bunches[:, 0]) / (binsize**2 * nshower)
#        if with_fov == 'n': # all photons
        h,edges = np.histogram(1e-2 * bunches[:, 2]
	    		   , bins = distances
	    		   , weights = weighted_pht
	    		   , range = [0., maxlen]
                           )
#        else: # 10 deg FoV
        wemis = np.sqrt(1 - bunches[:, 3]**2 - bunches[:, 4]**2)
        # Pointing telescope along the shower direction (on-axis) 
        # if theta(tel)!=theta_p -> off-axis observation
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180) 
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        weighted_pht_fov = weighted_pht[wemis >= fov]
        h_fov,edges = np.histogram(1e-2 * bunches_fov[:, 2]
				  , bins= distances
				  , weights= weighted_pht_fov
				  , range = [0., maxlen]
                                  )
    
    
    # Histogramming radially
    else:  
        r = np.sqrt((1e-2 * bunches[:, 1])**2 + (1e-2 * bunches[:, 2])**2)
        bunches = bunches[r < maxlen]
        r = r[r < maxlen]
        ring2 = ring[(r / (radius[1] - radius[0])).astype(int)]
        weighted_pht = np.abs(bunches[:, 0]) / (ring2 * nshower)
        # Store into an histogram
#        if with_fov == 'n':
        h, edges = np.histogram(r
				, bins = radius
				, weights = weighted_pht
				, range = [0., maxlen]
				)
#        else:
        wemis = np.sqrt(1 - bunches[:, 3]**2 - bunches[:, 4]**2)
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180) 
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        r_fov = r[wemis >= fov]
        ring2_fov = ring[(r_fov / (radius[1] - radius[0])).astype(int)]
        weighted_pht_fov = np.abs(bunches_fov[:, 0]) / (ring2_fov * nshower)
        h_fov, edges=np.histogram(r_fov
				  , bins = radius
				  , weights = weighted_pht_fov
				  , range = [0.,maxlen]
				  )
    return(h, h_fov)

# =========  DATA CARD =========
with open(sys.argv[2], 'r') as f:
    read_data = f.read()
    #print(read_data) 
    datacard = ascii.read(read_data,format='fixed_width_no_header'
			, delimiter=' '
                        ,col_starts=(0, 8, 39))
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

# Telescope pointing angle
onaxis = input("On-axis pointing (y/n)? ")
if onaxis == 'y':
    pointing_angle = theta_p
    pointing = 'onaxis'
else:
    pointing_angle = input("Off-axis angle? ")
    pointing = 'offaxis'

file = FortranFile(sys.argv[1], 'r')
# Histograms definition
binsize = 10 #meters
if x_area > y_area:
    maxlen  = 1e-2 * x_area / 2
    print("Histogramming along x-axis...")
    type_of_hist = 'x'
elif x_area < y_area:
    maxlen  = 1e-2 * y_area / 2
    print("Histogramming along y-axis...")
    type_of_hist = 'y'
else:
    type_of_hist = 'r'
    print("Histogramming radially...")
    maxlen  = 1e-2 * x_area / 2

numbins = int(maxlen / binsize)
breaks = numbins + 1
bunches = np.array([]).reshape(0, 7)
hist_c = np.zeros((2, numbins))
hist_f = np.zeros((2, numbins))
#with_fov='t'

if x_area != y_area:
    distances = np.linspace(0, maxlen, breaks)
    mids = ((distances[1] - distances[0]) / 2) + distances
    mids = mids[0:numbins]
else:
    radius = np.linspace(0, maxlen, breaks)
    mids = ((radius[1] - radius[0]) / 2) + radius
    mids = mids[0:numbins]
    ring = pi * ((radius[1:])**2 - (radius[0:numbins])**2)

# Control and debugging counters:
count = 0
lines = 0
photons = 0

while True:
    count = count + 1
    # Sort data in 21 subblocks of 39 lines each
    data = np.split(file.read_reals(dtype=np.float32).reshape(-1, 7), 21) 
    # It should be:
    # indices_boolean = [np.abs(i[0][0])<max(cersize,fluorsize) for i in data]
    indices_boolean = [np.abs(i[0][0]) < 100 for i in data] # select only subblocks of bunches
    indices = [i[0][0] for i in data] # store 1st element of each subblock
    bunches = np.vstack([bunches, np.vstack(compress(data, indices_boolean))])
    bunches = bunches[np.all(bunches != 0, axis = 1)] # drop those lines containing only zeros 

    if count == 10:
        lines = lines + len(bunches[bunches[:, 0] > 0]) + len(bunches[bunches[:,0] < 0])
        photons = photons + np.sum(bunches[:, 0])
        h_c = histogram(bunches[bunches[:, 0] > 0], x_area, y_area, pointing_angle)
        h_f = histogram(bunches[bunches[:, 0] < 0], x_area, y_area, pointing_angle)
        hist_c  = hist_c + h_c
        hist_f  = hist_f + h_f
        count = 0 # reset counter
        bunches = np.array([]).reshape(0, 7) # reset array
    
    if any( 3300 < i < 3303. for i in indices): #Flag indicating RUN END subblock
        if count<10:
            lines = lines + len(bunches[bunches[:, 0] > 0]) + len(bunches[bunches[:, 0] < 0])
            photons = photons + np.sum(bunches[:, 0])
            h_c = histogram(bunches[bunches[:, 0] > 0], x_area, y_area, pointing_angle)
            h_f = histogram(bunches[bunches[:, 0] < 0], x_area, y_area, pointing_angle)
            hist_c  = hist_c + h_c
            hist_f  = hist_f + h_f
            print(lines, ' lines (w/o zeros)')
            print(photons,' photons')
        file.close()
        break

np.savetxt('%iGeV_%ish_%ideg_%i%s_hist_%s.dat'%(E_prim,nshower,theta_p,pointing_angle,pointing,type_of_hist)
  ,np.transpose([mids, hist_c[0], hist_c[1], hist_f[0], hist_f[1]])
  ,newline = '\n'
  ,fmt = "%7.2f %1.6e %1.6e %1.6e %1.6e"
  ,header = (' Num_showers:%i \n E_primary (GeV): %i \n ID_prim_particle: %s \n Seeds: %i, %i \n'
              %(nshower,E_prim,prim_part,seed1,seed2) +
             ' Theta prim. part. incidence: %i deg \n Obs level (m): %i \n Atmosp model: %i'
              %(theta_p,obs_level ,atm_mod) +
             '\n Cerenk_bunch_size: %i \n Fluor_bunch_size: %i'
              %(cersize,fluorsize) +
             '\n  \n Distance to shower axis (m) | Phot_density_Cher/fluor (1/m2)'
             )
            )
print('Histogram stored into: %iGeV_%ish_%ideg_%i%s_hist_%s.dat'
       %(E_prim,nshower,theta_p,pointing_angle,pointing,type_of_hist)
      )


print("Reading completed")
print("-----------------","\n")
