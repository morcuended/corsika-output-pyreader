import numpy as np
from math import pi
from corsika.input import dcard


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
