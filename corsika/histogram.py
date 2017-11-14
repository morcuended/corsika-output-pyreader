

import numpy as np
from scipy.io import FortranFile
import pandas as pd
import matplotlib.pyplot as plt
from math import pi
from itertools import compress
from astropy.io import ascii
import matplotlib.style
import matplotlib as mpl
import sys, os


def PhotonBunches(bunches, x_area, y_area, theta, nshower):
    """
    Histogram photon bunches that reach observation level
    """
    # Histogram definition depending on the detection area
    binsize = 10  # meters

    if x_area > y_area:
        # print("Histogram along x-axis...")
        type_of_hist = 'x'
        maxlen = 1e-2 * x_area / 2
        numbins = int(maxlen / binsize)
        breaks = numbins + 1
        distances = np.linspace(0, maxlen, breaks)
        mids = ((distances[1] - distances[0]) / 2) + distances
        mids = mids[0:numbins]

    elif x_area < y_area:
        # print("Histogram along y-axis...")
        type_of_hist = 'y'
        maxlen = 1e-2 * y_area / 2
        numbins = int(maxlen / binsize)
        breaks = numbins + 1
        distances = np.linspace(0, maxlen, breaks)
        mids = ((distances[1] - distances[0]) / 2) + distances
        mids = mids[0:numbins]

    else:
        # print("Histogram radially...")
        type_of_hist = 'r'
        maxlen = 1e-2 * x_area / 2
        numbins = int(maxlen / binsize)
        breaks = numbins + 1
        radius = np.linspace(0, maxlen, breaks)
        mids = ((radius[1] - radius[0]) / 2) + radius
        mids = mids[0:numbins]
        ring = pi * ((radius[1:]) ** 2 - (radius[0:numbins]) ** 2)

    # bunches = np.array([]).reshape(0, 7)
    hist_c = np.zeros((2, numbins))
    hist_f = np.zeros((2, numbins))

    fov = np.cos(5 * pi / 180)  # FoV contraint (+/- 5 deg)
    # Histogramming along x-axis
    if x_area > y_area:
        weighted_pht = np.abs(bunches[:, 0]) / (binsize ** 2 * nshower)
        #        if with_fov == 'n': # all photons
        h, edges = np.histogram(1e-2 * bunches[:, 1],
                                bins=distances,
                                weights=weighted_pht,
                                range=[0., maxlen]
                                )
        wemis = np.sqrt(1 - bunches[:, 3] ** 2 - bunches[:, 4] ** 2)
        # Pointing the telescope along the shower direction (on-axis)
        # if theta(tel)!=theta_p -> off-axis observation
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180)
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        weighted_pht_fov = weighted_pht[wemis >= fov]
        h_fov, edges = np.histogram(1e-2 * bunches_fov[:, 1],
                                    bins = distances,
                                    weights = weighted_pht_fov,
                                    range = [0., maxlen]
                                    )
    elif x_area < y_area:
        weighted_pht = np.abs(bunches[:, 0]) / (binsize ** 2 * nshower)
        h, edges = np.histogram(1e-2 * bunches[:, 2],
                                bins = distances,
                                weights = weighted_pht,
                                range = [0., maxlen]
                                )
        wemis = np.sqrt(1 - bunches[:, 3] ** 2 - bunches[:, 4] ** 2)
        # Pointing telescope along the shower direction (on-axis)
        # if theta(tel)!=theta_p -> off-axis observation
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180)
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        weighted_pht_fov = weighted_pht[wemis >= fov]
        h_fov, edges = np.histogram(1e-2 * bunches_fov[:, 2],
                                    bins = distances,
                                    weights = weighted_pht_fov,
                                    range = [0., maxlen]
                                    )


    # Histogramming radially
    else:
        r = np.sqrt((1e-2 * bunches[:, 1]) ** 2 + (1e-2 * bunches[:, 2]) ** 2)
        bunches = bunches[r < maxlen]
        r = r[r < maxlen]
        ring2 = ring[(r / (radius[1] - radius[0])).astype(int)]
        weighted_pht = np.abs(bunches[:, 0]) / (ring2 * nshower)

        # Store into an histogram
        h, edges = np.histogram(r,
                                bins = radius,
                                weights = weighted_pht,
                                range = [0., maxlen]
                                )
        wemis = np.sqrt(1 - bunches[:, 3] ** 2 - bunches[:, 4] ** 2)
        wemis = wemis * np.cos(-theta * pi / 180) + bunches[:, 3] * np.sin(-theta * pi / 180)
        # where uemis=bunches[:,3]
        bunches_fov = bunches[wemis >= fov]
        r_fov = r[wemis >= fov]
        ring2_fov = ring[(r_fov / (radius[1] - radius[0])).astype(int)]
        weighted_pht_fov = np.abs(bunches_fov[:, 0]) / (ring2_fov * nshower)
        h_fov, edges = np.histogram(r_fov,
                                    bins = radius,
                                    weights = weighted_pht_fov,
                                    range = [0., maxlen]
                                    )
    return (h, h_fov)

