from scipy.io import FortranFile
from itertools import compress
import argparse
import os
from corsika.input import dcard
from corsika.histogram import hist, hist2d, type_of_hist, hist_c, hist_f, hist_2d, mids
import numpy as np


"""
Histogram Cherenkov and fluorescence photon bunches on the ground.
In both ways, 2D histogram or lateral photon density profile.
"""


parser = argparse.ArgumentParser(description='Histogram photon density '
                                             'profile from CERnnnnnn '
                                             'standard binary CORSIKA '
                                             'output.',
                                 usage='reader.py path/to/CORSIKA/output/'
                                       'files [pointing] [hist2d]')

parser.add_argument('Directory',
                    help='path to directory containing output CORSIKA data',
                    nargs='?', default=os.getcwd())

parser.add_argument('-p', '--pointing',
                    type=float,
                    help='Telescope pointing angle with respect '
                         'to vertical direction',
                    default=0.0)

parser.add_argument('-f', '--fov',
                    type=float,
                    help='Diameter field of view (FoV) of the telescope',
                    default=10.0)

parser.add_argument('-h2d', '--hist2d',
                    help='Generate 2D histogram of photon density '
                         ' on the ground',
                    action="store_true")

args = parser.parse_args()


print('\n-------  CORSIKA reader  -------')

os.chdir(args.Directory)

print('You are in: ', os.getcwd())


# Telescope pointing
if args.pointing == dcard['THETAP']:
    pointing_angle = dcard['THETAP']
    pointing = 'onaxis'
    off_axis = 0.
else:
    pointing_angle = args.pointing
    pointing = 'offaxis'
    off_axis = pointing_angle - dcard['THETAP']

# Open CORSIKA output binary file
file = FortranFile("CER000001", 'r')


# Control counter and empty photon bunches array definition:
count = 0
bunches = np.array([]).reshape(0, 7)

while True:
    count = count + 1
    # Sort data in 21 sub-blocks of 39 lines each
    data = np.split(file.read_reals(dtype=np.float32).reshape(-1, 7), 21)
    # Select only sub-blocks of bunches
    indices_boolean = [np.abs(i[0][0]) < 3000 for i in data]
    # Store 1st element of each sub-block
    indices = [i[0][0] for i in data]
    bunches = np.vstack([bunches, np.vstack(compress(data, indices_boolean))])
    # Drop those lines containing only zeros
    bunches = bunches[np.all(bunches != 0, axis=1)]

    if count == 10:
        if not args.hist2d:
            h_c = hist(bunches[bunches[:, 0] > 0],
                       dcard['XCERARY'],
                       dcard['YCERARY'],
                       pointing_angle,
                       dcard['NSHOW'],
                       args.fov
                       )
            h_f = hist(bunches[bunches[:, 0] < 0],
                       dcard['XCERARY'],
                       dcard['YCERARY'],
                       pointing_angle,
                       dcard['NSHOW'],
                       args.fov
                       )
            hist_c = hist_c + h_c
            hist_f = hist_f + h_f
        if args.hist2d:
            h2d = hist2d(bunches,
                         dcard['NSHOW']
                         )
            hist_2d = hist_2d + h2d
        count = 0  # reset counter
        bunches = np.array([]).reshape(0, 7)  # reset array

    if any(3300 < i < 3303. for i in indices):  # RUN END subblock
        if count < 10:
            if not args.hist2d:
                h_c = hist(bunches[bunches[:, 0] > 0],
                           dcard['XCERARY'],
                           dcard['YCERARY'],
                           pointing_angle,
                           dcard['NSHOW'],
                           args.fov
                           )
                h_f = hist(bunches[bunches[:, 0] < 0],
                           dcard['XCERARY'],
                           dcard['YCERARY'],
                           pointing_angle,
                           dcard['NSHOW'],
                           args.fov
                           )
                hist_c = hist_c + h_c
                hist_f = hist_f + h_f
            if args.hist2d:
                h2d = hist2d(bunches,
                             dcard['NSHOW']
                             )
                hist_2d = hist_2d + h2d
        file.close()
        break

if args.hist2d:
    np.savetxt('%iGeV_%ish_%ideg_hist2d.dat' % (dcard['ERANGE'],
                                                dcard['NSHOW'],
                                                dcard['THETAP']),
               hist_2d,
               fmt="%1.7e"
               )
    print('2D Histogram stored into: %iGeV_%ish_%ideg_hist2d.dat' %
          (dcard['ERANGE'], dcard['NSHOW'], dcard['THETAP'])
          )

if not args.hist2d:
    np.savetxt('%iGeV_%ish_%ideg_%i%s_hist_%s_%iFoV.dat' % (dcard['ERANGE'],
                                                            dcard['NSHOW'],
                                                            dcard['THETAP'],
                                                            abs(off_axis),
                                                            pointing,
                                                            type_of_hist,
                                                            args.fov),
               np.transpose([mids, hist_c[0], hist_c[1],  # Cherenkov
                             hist_f[0], hist_f[1]]),  # fluorescence
               newline='\n',
               fmt="%7.2f %1.6e %1.6e %1.6e %1.6e",
               header=(' Number of showers: %i \n Energy primary [GeV]: %i \n'
                       ' ID primary particle: %s \n'
                       ' Seeds for sequences 1, 2, 3: %i, %i, %i \n'
                       % (dcard['NSHOW'],
                          dcard['ERANGE'],
                          dcard['PRMPAR'],
                          dcard['SEED'], dcard['SEED1'], dcard['SEED2'])
                       +
                       ' Theta primary particle [deg]: %i \n'
                       ' Observation level [cm]: %i \n'
                       ' Atm. model: %i \n'
                       ' Diameter FoV [deg]: %i \n'
                       % (dcard['THETAP'],
                          dcard['OBSLEV'],
                          dcard['ATMOD'],
                          args.fov)
                       +
                       ' Cherenkov bunch size: %i \n Fluor bunch size: %i \n'
                       % (dcard['CERSIZ'], dcard['FLSIZE'])
                       +
                       ' Distance to core [m] |'
                       ' Ph. density [1/m2] C (all) | C (FoV) |'
                       ' F (all) | F (FoV)'
                       )
               )
    print('Histogram stored into: %iGeV_%ish_%ideg_%i%s_hist_%s_%iFoV.dat \n' %
          (dcard['ERANGE'], dcard['NSHOW'], dcard['THETAP'],
           abs(off_axis), pointing, type_of_hist, args.fov) +
          'Telescope pointing angle [deg]: %4.1f \n' % pointing_angle +
          'Off-axis angle wrt the shower core [deg]: %4.1f \n' % off_axis +
          'Diameter FoV [deg]: %3.1f' % args.fov
          )

print("--------------------------------\n")
