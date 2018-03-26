from astropy.io import ascii


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
