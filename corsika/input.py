from astropy.io import ascii
import re


#  Load Input CORSIKA Card

dcard = {}
times = 0

with open("all", 'r') as f:
    """
    Read input variables from DATACARD
    """
    read_data = f.read()
    raw_datacard = ascii.read(read_data, format='fixed_width_no_header',
                              delimiter=' ', col_starts=(0, 8, 37))

for keyword, value, description in raw_datacard:
    while True:
        try:
            if not str(keyword) in dcard:
                dcard[str(keyword)] = [float(x) for x in value.split()]
            else:
                times += 1
                dcard[str(keyword)+str(times)] = [float(x) for x in value.split()]
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
dcard['SEED1'] = dcard['SEED1'][0]
dcard['SEED2'] = dcard['SEED2'][0]
dcard['NSHOW'] = dcard['NSHOW'][0]
dcard['OBSLEV'] = dcard['OBSLEV'][0]
dcard['CERSIZ'] = dcard['CERSIZ'][0]
dcard['FLSIZE'] = dcard['FLSIZE'][0]
dcard['PRMPAR'] = dcard['PRMPAR'][0]
dcard['ATMOD'] = dcard['ATMOD'][0]


#  Load DATnnnnnn.dbase file

try:
    filename = "DAT000001.dbase"

    with open(filename) as f:
            content = f.readlines()

    RunHeader = {}
    for line in content:
        raw_dbase = re.search('#(.*)', line)
        dbase = raw_dbase.group(1).split("#")
        for keyword, value in zip(dbase[0::2], dbase[1::2]):
            try:
                RunHeader[str(keyword)] = float(value)
            except ValueError:
                RunHeader[str(keyword)] = str(value)
except FileNotFoundError:
    print("DATnnnnnn.dbase file not found")
