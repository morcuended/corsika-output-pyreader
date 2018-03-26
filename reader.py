from scipy.io import FortranFile
import numpy as np
from math import pi
from itertools import compress
from astropy.io import ascii
import argparse
import os


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
                    help='path to directory containg output CORSIKA data',
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


def histogram(photon_bunches, x_area, y_area, pointing_ang, nshower):
    """
    :param photon_bunches: array containing only photon bunches
    sub-blocks. The information contained is:
    [N_phot, x, y, uemis, vemis, t_emis, height of emission]
    :param x_area: x dimension of detection area (in cm)
    :param y_area: x dimension of detection area (in cm)
    :param pointing_ang: telescope pointing angle
                         with respect to vertical direction
    :param nshower: number of simulated showers
    :return: photons histogrammed radially/along x or y-axis
    """
    fov = np.cos(0.5 * args.fov * pi / 180)  # FoV constraint radially
    # Histogram along x-axis
    if x_area > y_area:
        weighted_pht = np.abs(photon_bunches[:, 0]) / (binsize**2 * nshower)
        h, edges = np.histogram(1e-2 * photon_bunches[:, 1],
                                bins=distances,
                                weights=weighted_pht,
                                range=[0., maxlen]
                                )
        wemis = np.sqrt(1-photon_bunches[:, 3]**2 - photon_bunches[:, 4]**2)
        # Pointing telescope along the pointing angle direction
        # set as input argument
        wemis = (wemis * np.cos(-pointing_ang * pi / 180)
                 + photon_bunches[:, 3] * np.sin(-pointing_ang * pi / 180))
        bunches_fov = photon_bunches[wemis >= fov]
        weighted_pht_fov = weighted_pht[wemis >= fov]
        h_fov, edges = np.histogram(1e-2 * bunches_fov[:, 1],
                                    bins=distances,
                                    weights=weighted_pht_fov,
                                    range=[0., maxlen]
                                    )
    elif x_area < y_area:
        weighted_pht = np.abs(photon_bunches[:, 0]) / (binsize**2 * nshower)
        h, edges = np.histogram(1e-2 * photon_bunches[:, 2],
                                bins=distances,
                                weights=weighted_pht,
                                range=[0., maxlen]
                                )
        wemis = np.sqrt(1 - photon_bunches[:, 3]**2 - photon_bunches[:, 4]**2)
        wemis = (wemis * np.cos(-pointing_ang * pi / 180)
                 + photon_bunches[:, 3] * np.sin(-pointing_ang * pi / 180))
        bunches_fov = photon_bunches[wemis >= fov]
        weighted_pht_fov = weighted_pht[wemis >= fov]
        h_fov, edges = np.histogram(1e-2 * bunches_fov[:, 2],
                                    bins=distances,
                                    weights=weighted_pht_fov,
                                    range=[0., maxlen]
                                    )

    # Histogram radially
    else:
        r = np.sqrt((1e-2 * photon_bunches[:, 1])**2
                    + (1e-2 * photon_bunches[:, 2])**2)
        photon_bunches = photon_bunches[r < maxlen]
        r = r[r < maxlen]
        ring2 = ring[(r / (radius[1] - radius[0])).astype(int)]
        weighted_pht = np.abs(photon_bunches[:, 0]) / (ring2 * nshower)

        # Store into the histogram
        h, edges = np.histogram(r,
                                bins=radius,
                                weights=weighted_pht,
                                range=[0., maxlen]
                                )
        wemis = np.sqrt(1 - photon_bunches[:, 3]**2 - photon_bunches[:, 4]**2)
        wemis = (wemis * np.cos(-pointing_ang * pi / 180)
                 + photon_bunches[:, 3] * np.sin(-pointing_ang * pi / 180))
        bunches_fov = photon_bunches[wemis >= fov]
        r_fov = r[wemis >= fov]
        ring2_fov = ring[(r_fov / (radius[1] - radius[0])).astype(int)]
        weighted_pht_fov = np.abs(bunches_fov[:, 0]) / (ring2_fov * nshower)
        h_fov, edges = np.histogram(r_fov,
                                    bins=radius,
                                    weights=weighted_pht_fov,
                                    range=[0., maxlen]
                                    )
    return h, h_fov


def hist2d(photon_bunches, nshower):
    """
    :param photon_bunches: array containing only photon bunches
    sub-blocks. The information contained is:
    [N_phot, x, y, uemis, vemis, t_emis, height of emission]
    :param nshower: number of simulated showers
    :return: photons  radially/along x or y-axis
    """
    # 2D Histogram
    weighted_pht = np.abs(photon_bunches[:, 0]) / (binsize**2 * nshower)
    histogram2d, xedges, yedges = np.histogram2d(1e-2 * photon_bunches[:, 1],
                                                 1e-2 * photon_bunches[:, 2],
                                                 bins=distances2D,
                                                 weights=weighted_pht,
                                                 range=([-maxlen, maxlen],
                                                        [-maxlen, maxlen])
                                                 )
    return histogram2d


dcard = {}
with open("all", 'r') as f:
    """
    Read input variables from DATACARD
    """
    read_data = f.read()
    raw_datacard = ascii.read(read_data, format='fixed_width_no_header',
                              delimiter=' ', col_starts=(0, 8, 39))

for keyword, value, description in raw_datacard:
    while True:
        try:
            dcard[str(keyword)] = [float(x) for x in value.split()]
        except ValueError:
            dcard[str(keyword)] = [str(x) for x in value.split()]
        except AttributeError:
            print('DataCard successfully loaded')
        break

# Redefining some of the fields
dcard['THETAP'] = dcard['THETAP'][0]
dcard['XCERARY'] = dcard['CERARY'][4]
dcard['YCERARY'] = dcard['CERARY'][5]
dcard['ERANGE'] = dcard['ERANGE'][0]
dcard['SEED'] = dcard['SEED'][0]
dcard['NSHOW'] = dcard['NSHOW'][0]
dcard['OBSLEV'] = dcard['OBSLEV'][0]
dcard['CERSIZ'] = dcard['CERSIZ'][0]
dcard['FLSIZE'] = dcard['FLSIZE'][0]
dcard['PRMPAR'] = dcard['PRMPAR'][0]
dcard['ATMOD'] = dcard['ATMOD'][0]


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

# Histograms definition
binsize = 10  # meters
maxlen = 0.5e-2 * max(dcard['XCERARY'], dcard['YCERARY'])
numbins = int(maxlen / binsize)
breaks = numbins + 1
distances = np.linspace(0, maxlen, breaks)
distances2D = np.linspace(-maxlen, maxlen, 2*breaks-1)
mids = ((distances[1] - distances[0]) / 2) + distances
mids = mids[0:numbins]
hist_c = np.zeros((2, numbins))
hist_f = np.zeros((2, numbins))
Histogram2D = np.zeros((2*numbins, 2*numbins))


if dcard['XCERARY'] > dcard['YCERARY']:
    type_of_hist = 'x'

elif dcard['XCERARY'] < dcard['YCERARY']:
    type_of_hist = 'y'

else:
    type_of_hist = 'r'
    radius = np.linspace(0, maxlen, breaks)
    mids = ((radius[1] - radius[0]) / 2) + radius
    mids = mids[0:numbins]
    ring = pi * ((radius[1:]) ** 2 - (radius[0:numbins]) ** 2)

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
            h_c = histogram(bunches[bunches[:, 0] > 0],
                            dcard['XCERARY'],
                            dcard['YCERARY'],
                            pointing_angle,
                            dcard['NSHOW']
                            )
            h_f = histogram(bunches[bunches[:, 0] < 0],
                            dcard['XCERARY'],
                            dcard['YCERARY'],
                            pointing_angle,
                            dcard['NSHOW']
                            )
            hist_c = hist_c + h_c
            hist_f = hist_f + h_f
        if args.hist2d:
            h2d = hist2d(bunches,
                         dcard['NSHOW']
                         )
            Histogram2D = Histogram2D + h2d
        count = 0  # reset counter
        bunches = np.array([]).reshape(0, 7)  # reset array

    if any(3300 < i < 3303. for i in indices):  # RUN END subblock
        if count < 10:
            if not args.hist2d:
                h_c = histogram(bunches[bunches[:, 0] > 0],
                                dcard['XCERARY'],
                                dcard['YCERARY'],
                                pointing_angle,
                                dcard['NSHOW']
                                )
                h_f = histogram(bunches[bunches[:, 0] < 0],
                                dcard['XCERARY'],
                                dcard['YCERARY'],
                                pointing_angle,
                                dcard['NSHOW']
                                )
                hist_c = hist_c + h_c
                hist_f = hist_f + h_f
            if args.hist2d:
                h2d = hist2d(bunches,
                             dcard['NSHOW']
                             )
                Histogram2D = Histogram2D + h2d
        file.close()
        break

if args.hist2d:
    np.savetxt('%iGeV_%ish_%ideg_hist2d.dat' % (dcard['ERANGE'],
                                                dcard['NSHOW'],
                                                dcard['THETAP']),
               Histogram2D,
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
                       ' Seed for Cherenkov emission: %i \n'
                       % (dcard['NSHOW'],
                          dcard['ERANGE'],
                          dcard['PRMPAR'],
                          dcard['SEED'])
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
