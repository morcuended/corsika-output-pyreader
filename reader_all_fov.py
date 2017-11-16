from scipy.io import FortranFile
import numpy as np
from math import pi
from itertools import compress
from astropy.io import ascii
import sys

''' 
Histogram Cherenkov and fluorescence photon bunches
on the ground.
 
Usage:

    python reader.py path_to_output_file_from_CORSIKA DATACARD_file

TO DO list:
  (x) Establish a way of pointing the telescope to any direction
      via theta angle. 
  (x) Ask if user wants on-axis pointing, if not, ask for the 
      pointing angle of the telescope.
  (-) Convert pointing angle to off-axis angle with respect to the 
      shower axis direction.
  (-) If any argument is passed, CER000001 should be set
      as default input file.
  (-) Average x/y-histogram.The idea is to have to different DATACARD
      with detector area defined along both axis (x/y) with the
      corresponding shower direction for each one.
  (-) Check which dimension is the largest one, then histogram along
      that direction.
'''

print('\n-------  CORSIKA reader  -------')


def histogram(bunches, x_area, y_area, theta, nshower):
    fov = np.cos(5 * pi / 180)  # FoV constraint (+/- 5 deg)
    # Histogram along x-axis
    if x_area > y_area:  
        weighted_pht = np.abs(bunches[:, 0]) / (binsize**2 * nshower)
        h, edges = np.histogram(1e-2 * bunches[:, 1],
                                bins=distances,
                                weights=weighted_pht,
                                range=[0., maxlen]
                                )
#        else: # 10 deg FoV
        wemis = np.sqrt(1-bunches[:, 3]**2 - bunches[:, 4]**2)
        # Pointing the telescope along the shower direction (on-axis) 
        # if theta(tel)!=theta_p -> off-axis observation
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180) 
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        weighted_pht_fov = weighted_pht[wemis >= fov]
        h_fov, edges = np.histogram(1e-2 * bunches_fov[:, 1],
                                    bins=distances,
                                    weights=weighted_pht_fov,
                                    range=[0., maxlen]
                                    )
    elif x_area < y_area:  
        weighted_pht = np.abs(bunches[:, 0]) / (binsize**2 * nshower)
        h, edges = np.histogram(1e-2 * bunches[:, 2],
                                bins=distances,
                                weights=weighted_pht,
                                range=[0., maxlen]
                                )
        wemis = np.sqrt(1 - bunches[:, 3]**2 - bunches[:, 4]**2)
        # Pointing telescope along the shower direction (on-axis) 
        # if theta(tel)!=theta_p -> off-axis observation
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180) 
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        weighted_pht_fov = weighted_pht[wemis >= fov]
        h_fov, edges = np.histogram(1e-2 * bunches_fov[:, 2],
                                    bins=distances,
                                    weights=weighted_pht_fov,
                                    range=[0., maxlen]
                                    )
    
    # Histogram radially
    else:  
        r = np.sqrt((1e-2 * bunches[:, 1])**2 + (1e-2 * bunches[:, 2])**2)
        bunches = bunches[r < maxlen]
        r = r[r < maxlen]
        ring2 = ring[(r / (radius[1] - radius[0])).astype(int)]
        weighted_pht = np.abs(bunches[:, 0]) / (ring2 * nshower)
        # Store into an histogram
        h, edges = np.histogram(r,
                                bins=radius,
                                weights=weighted_pht,
                                range=[0., maxlen]
                                )
        wemis = np.sqrt(1 - bunches[:, 3]**2 - bunches[:, 4]**2)
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180) 
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        r_fov = r[wemis >= fov]
        ring2_fov = ring[(r_fov / (radius[1] - radius[0])).astype(int)]
        weighted_pht_fov = np.abs(bunches_fov[:, 0]) / (ring2_fov * nshower)
        h_fov, edges = np.histogram(r_fov,
                                    bins=radius,
                                    weights=weighted_pht_fov,
                                    range=[0., maxlen]
                                    )
    return h, h_fov


with open(sys.argv[2], 'r') as f:
    """
    Read input variables from DATACARD
    """
    read_data = f.read()
    datacard = ascii.read(read_data, format='fixed_width_no_header',
                          delimiter=' ', col_starts=(0, 8, 39))

data_card = dict(NSHOW=int(datacard[2][1].split()[0]),
                 ERANGE=float(datacard[5][1].split()[1]),
                 PRMPAR=int(datacard[3][1].split()[0]),
                 SEED1=int(datacard[8][1].split()[0]),
                 SEED2=int(datacard[9][1].split()[0]),
                 THETAP=float(datacard[6][1].split()[0]),
                 PHIP=float(datacard[7][1].split()[0]),
                 OBSLEV=float(datacard[10][1].split()[0]) * 1e-2,
                 ATMOD=int(datacard[18][1].split()[0]),
                 CERSIZ=float(datacard[23][1].split()[0]),
                 FLSIZE=float(datacard[24][1].split()[0]),
                 XCERARY=float(datacard[27][1].split()[4]),
                 YCERARY=float(datacard[27][1].split()[5])
                 )

# Telescope pointing angle
onaxis = input("On-axis pointing (y/n)? ")
if onaxis == 'y':
    pointing_angle = data_card['THETAP']
    pointing = 'onaxis'
else:
    pointing_angle = input("Off-axis angle? ")
    pointing = 'offaxis'

# Open CORSIKA output binary file
file = FortranFile(sys.argv[1], 'r')

# Histograms definition
binsize = 10  # meters

if data_card['XCERARY'] > data_card['YCERARY']:
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
    maxlen = 1e-2 * data_card['XCERARY'] / 2
    numbins = int(maxlen / binsize)
    breaks = numbins + 1
    radius = np.linspace(0, maxlen, breaks)
    mids = ((radius[1] - radius[0]) / 2) + radius
    mids = mids[0:numbins]
    ring = pi * ((radius[1:]) ** 2 - (radius[0:numbins]) ** 2)

bunches = np.array([]).reshape(0, 7)
hist_c = np.zeros((2, numbins))
hist_f = np.zeros((2, numbins))

# Control and debugging counters:
count = 0


while True:
    count = count + 1
    # Sort data in 21 sub-blocks of 39 lines each
    data = np.split(file.read_reals(dtype=np.float32).reshape(-1, 7), 21) 
    # It should be:
    # indices_boolean = [np.abs(i[0][0])<max(cersize, fluorsize) for i in data]
    # Select only sub-blocks of bunches
    indices_boolean = [np.abs(i[0][0]) < 100 for i in data]
    # Store 1st element of each sub-block
    indices = [i[0][0] for i in data]
    bunches = np.vstack([bunches, np.vstack(compress(data, indices_boolean))])
    # Drop those lines containing only zeros
    bunches = bunches[np.all(bunches != 0, axis=1)]

    if count == 10:
        h_c = histogram(bunches[bunches[:, 0] > 0],
                        data_card['XCERARY'],
                        data_card['YCERARY'],
                        pointing_angle,
                        data_card['NSHOW']
                        )
        h_f = histogram(bunches[bunches[:, 0] < 0],
                        data_card['XCERARY'],
                        data_card['YCERARY'],
                        pointing_angle,
                        data_card['NSHOW']
                        )
        hist_c = hist_c + h_c
        hist_f = hist_f + h_f
        count = 0  # reset counter
        bunches = np.array([]).reshape(0, 7)  # reset array
    
    if any(3300 < i < 3303. for i in indices):  # Flag indicating RUN END subblock
        if count < 10:
            h_c = histogram(bunches[bunches[:, 0] > 0],
                            data_card['XCERARY'],
                            data_card['YCERARY'],
                            pointing_angle,
                            data_card['NSHOW']
                            )
            h_f = histogram(bunches[bunches[:, 0] < 0],
                            data_card['XCERARY'],
                            data_card['YCERARY'],
                            pointing_angle,
                            data_card['NSHOW']
                            )
            hist_c = hist_c + h_c
            hist_f = hist_f + h_f
        file.close()
        break

np.savetxt('%iGeV_%ish_%ideg_%i%s_hist_%s.dat' % (data_card['ERANGE'],
                                                  data_card['NSHOW'],
                                                  data_card['THETAP'],
                                                  pointing_angle,
                                                  pointing,
                                                  type_of_hist),
           np.transpose([mids, hist_c[0], hist_c[1], hist_f[0], hist_f[1]]),
           newline='\n',
           fmt="%7.2f %1.6e %1.6e %1.6e %1.6e",
           header=(' Num_showers:%i \n E_primary (GeV): %i \n ID_prim_particle: %s \n Seeds: %i, %i \n'
                   % (data_card['NSHOW'],
                      data_card['ERANGE'],
                      data_card['PRMPAR'],
                      data_card['SEED1'],
                      data_card['SEED2'])
                   +
                   ' Theta prim. part. incidence: %i deg \n Obs level (m): %i \n Atmosp model: %i'
                   % (data_card['THETAP'],
                      data_card['OBSLEV'],
                      data_card['ATMOD'])
                   +
                   '\n Cerenk_bunch_size: %i \n Fluor_bunch_size: %i'
                   % (data_card['CERSIZ'], data_card['FLSIZE'])
                   +
                   '\n Distance to shower axis (m) | Phot_density_Cher/fluor (1/m2)'
                   )
           )
print('Histogram stored into: %iGeV_%ish_%ideg_%i%s_hist_%s.dat' %
      (data_card['ERANGE'], data_card['NSHOW'], data_card['THETAP'],
       pointing_angle, pointing, type_of_hist)
      )

print("--------------------------------\n")
