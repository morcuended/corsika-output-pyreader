from corsika.histogram_definition import *
import numpy as np
from math import pi


""""
Histograms definition
"""

bin_size = 10  # meters
max_len = 0.5e-2 * max(dcard['XCERARY'], dcard['YCERARY'])
num_bins = int(max_len / bin_size)
breaks = num_bins + 1
distances = np.linspace(0, max_len, breaks)
distances2d = np.linspace(-max_len, max_len, 2 * breaks - 1)
mids = ((distances[1] - distances[0]) / 2) + distances
mids = mids[0:num_bins]

if dcard['XCERARY'] > dcard['YCERARY']:
    type_of_hist = 'x'

elif dcard['XCERARY'] < dcard['YCERARY']:
    type_of_hist = 'y'

else:
    type_of_hist = 'r'
    radius = np.linspace(0, max_len, breaks)
    mids = ((radius[1] - radius[0]) / 2) + radius
    mids = mids[0:num_bins]
    ring = pi * ((radius[1:]) ** 2 - (radius[0:num_bins]) ** 2)

hist_c = np.zeros((2, num_bins))
hist_f = np.zeros((2, num_bins))
hist_2d = np.zeros((2 * num_bins, 2 * num_bins))


def hist(photon_bunches, x_detector, y_detector, pointing_ang, n_shower, fov_diameter):
    """
    :param photon_bunches: array containing only photon bunches
    sub-blocks. The information contained is:
    [N_phot, x, y, u_emis, v_emis, t_emis, height of emission]
    :param x_detector: x dimension of detection area (in cm)
    :param y_detector: x dimension of detection area (in cm)
    :param pointing_ang: telescope pointing angle
                         with respect to vertical direction
    :param n_shower: number of simulated showers
    :param fov_diameter: fov (angular diameter) of telescope
    :return: photons histogrammed radially/along x or y-axis
    """

    fov = np.cos(0.5 * fov_diameter * pi / 180)  # FoV constraint radially

    # Histogram along x-axis
    if x_detector > y_detector:
        weighted_pht = np.abs(photon_bunches[:, 0]) / (bin_size ** 2 * n_shower)
        h, edges = np.histogram(1e-2 * photon_bunches[:, 1],
                                bins=distances,
                                weights=weighted_pht,
                                range=[0., max_len]
                                )
        w_emis = np.sqrt(1 - photon_bunches[:, 3] ** 2 - photon_bunches[:, 4] ** 2)
        # Pointing telescope along the pointing angle direction
        # set as input argument
        w_emis = (w_emis * np.cos(-pointing_ang * pi / 180)
                  + photon_bunches[:, 3] * np.sin(-pointing_ang * pi / 180))
        bunches_fov = photon_bunches[w_emis >= fov]
        weighted_pht_fov = weighted_pht[w_emis >= fov]
        h_fov, edges = np.histogram(1e-2 * bunches_fov[:, 1],
                                    bins=distances,
                                    weights=weighted_pht_fov,
                                    range=[0., max_len]
                                    )
    elif x_detector < y_detector:
        weighted_pht = np.abs(photon_bunches[:, 0]) / (bin_size ** 2 * n_shower)
        h, edges = np.histogram(1e-2 * photon_bunches[:, 2],
                                bins=distances,
                                weights=weighted_pht,
                                range=[0., max_len]
                                )
        w_emis = np.sqrt(1 - photon_bunches[:, 3] ** 2 - photon_bunches[:, 4] ** 2)
        w_emis = (w_emis * np.cos(-pointing_ang * pi / 180)
                  + photon_bunches[:, 3] * np.sin(-pointing_ang * pi / 180))
        bunches_fov = photon_bunches[w_emis >= fov]
        weighted_pht_fov = weighted_pht[w_emis >= fov]
        h_fov, edges = np.histogram(1e-2 * bunches_fov[:, 2],
                                    bins=distances,
                                    weights=weighted_pht_fov,
                                    range=[0., max_len]
                                    )

    # Histogram radially
    else:
        r = np.sqrt((1e-2 * photon_bunches[:, 1]) ** 2
                    + (1e-2 * photon_bunches[:, 2]) ** 2)
        photon_bunches = photon_bunches[r < max_len]
        r = r[r < max_len]
        ring2 = ring[(r / (radius[1] - radius[0])).astype(int)]
        weighted_pht = np.abs(photon_bunches[:, 0]) / (ring2 * n_shower)

        # Store into the histogram
        h, edges = np.histogram(r,
                                bins=radius,
                                weights=weighted_pht,
                                range=[0., max_len]
                                )
        w_emis = np.sqrt(1 - photon_bunches[:, 3] ** 2 - photon_bunches[:, 4] ** 2)
        w_emis = (w_emis * np.cos(-pointing_ang * pi / 180)
                  + photon_bunches[:, 3] * np.sin(-pointing_ang * pi / 180))
        bunches_fov = photon_bunches[w_emis >= fov]
        r_fov = r[w_emis >= fov]
        ring2_fov = ring[(r_fov / (radius[1] - radius[0])).astype(int)]
        weighted_pht_fov = np.abs(bunches_fov[:, 0]) / (ring2_fov * n_shower)
        h_fov, edges = np.histogram(r_fov,
                                    bins=radius,
                                    weights=weighted_pht_fov,
                                    range=[0., max_len]
                                    )
    return h, h_fov


def hist2d(photon_bunches, n_shower):
    """
    :param photon_bunches: array containing only photon bunches
    sub-blocks. The information contained is:
    [N_phot, x, y, u_emis, v_emis, t_emis, height of emission]
    :param n_shower: number of simulated showers
    :return: photons  radially/along x or y-axis
    """
    # 2D Histogram
    weighted_pht = np.abs(photon_bunches[:, 0]) / (bin_size ** 2 * n_shower)
    histogram2d, xedges, yedges = np.histogram2d(1e-2 * photon_bunches[:, 1],
                                                 1e-2 * photon_bunches[:, 2],
                                                 bins=distances2d,
                                                 weights=weighted_pht,
                                                 range=([-max_len, max_len],
                                                        [-max_len, max_len])
                                                 )
    return histogram2d
