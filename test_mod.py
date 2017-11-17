"""
This script reads Cherenkov and fluorescence photon bunches
from all events in the CORSIKA output.

Usage: It takes CERnnnnnn and DATACARD files as input arguments:

    python test.py CERnnnnnn all

Variables in capital letters corespond to datacard input arguments
"""

import histogram_mod
import DataCard
import sys
from itertools import compress
import numpy as np
from scipy.io import FortranFile

data_card = DataCard.read(sys.argv[2])

# Define telescope pointing angle
onaxis = input("On-axis pointing (y/n)? ")
if onaxis == 'y':
    pointing_angle = data_card['THETAP']
    pointing = 'onaxis'
else:
    pointing_angle = input("Off-axis angle? ")
    pointing = 'offaxis'

# Definition of the histogram depending on the detection area
bunches = np.array([]).reshape(0, 7)
binsize = 10  # meters
type_of_hist, hist_c, hist_f, mids = histogram_mod.PhotonBunches(bunches[bunches[:, 0] > 0],
                                                                 data_card['XCERARY'],
                                                                 data_card['YCERARY'],
                                                                 pointing_angle,
                                                                 data_card['NSHOW'],
                                                                 definition=True
                                                                 )

# Open binary CORSIKA output file
file = FortranFile(sys.argv[1], 'r')

# Control and debugging counters:
count = 0


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
        h_c = histogram_mod.PhotonBunches(bunches[bunches[:, 0] > 0],
                                          data_card['XCERARY'],
                                          data_card['YCERARY'],
                                          pointing_angle,
                                          data_card['NSHOW']
                                          )
        h_f = histogram_mod.PhotonBunches(bunches[bunches[:, 0] < 0],
                                          data_card['XCERARY'],
                                          data_card['YCERARY'],
                                          pointing_angle,
                                          data_card['NSHOW']
                                          )
        hist_c = hist_c + h_c
        hist_f = hist_f + h_f
        count = 0  # reset counter
        bunches = np.array([]).reshape(0, 7)  # reset array

    if any(3300 < i < 3303. for i in indices):  # Flag indicating RUN END sub-block
        if count < 10:
            h_c = histogram_mod.PhotonBunches(bunches[bunches[:, 0] > 0],
                                              data_card['XCERARY'],
                                              data_card['YCERARY'],
                                              pointing_angle,
                                              data_card['NSHOW']
                                              )
            h_f = histogram_mod.PhotonBunches(bunches[bunches[:, 0] < 0],
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
