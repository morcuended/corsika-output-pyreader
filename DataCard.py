from astropy.io import ascii

import os
import sys

# Open DATA CARD
def read(data_card_file):
    """
    It reads the input arguments passed in the DATA CARD.

    So far not implemented for all of the input arguments.

    In this configuration, "all" file must be placed in the
    parent directory, from where test.py script is run.

    It returns only one dictionary storing all DATA CARD variables
    """
    with open(data_card_file, 'r') as f:
        read_data = f.read()
        datacard = ascii.read(read_data, format='fixed_width_no_header',
                              delimiter=' ', col_starts=(0, 8, 39)
                              )
    data_card = {
            'NSHOW': int(datacard[2][1].split()[0]),
            'ERANGE': float(datacard[5][1].split()[1]),
            'PRMPAR': int(datacard[3][1].split()[0]),
            'SEED1': int(datacard[8][1].split()[0]),
            'SEED2': int(datacard[9][1].split()[0]),
            'THETAP': float(datacard[6][1].split()[0]),
            'PHIP': float(datacard[7][1].split()[0]),
            'OBSLEV': float(datacard[10][1].split()[0]) * 1e-2,
            'ATMOD': int(datacard[18][1].split()[0]),
            'CERSIZ': float(datacard[23][1].split()[0]),
            'FLSIZE': float(datacard[24][1].split()[0]),
            'XCERARY': float(datacard[27][1].split()[4]),
            'YCERARY': float(datacard[27][1].split()[5])
    }

    return data_card
