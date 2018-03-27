# CORSIKA output reader

Python-based script for reading and histogramming Cherenkov/Fluorescence photons from 
CORSIKA standard output *CERnnnnnn* (without IACT option). It also performs 2D histograms of both light components on the ground.

Tested for [CORSIKA Version 7.6300](https://web.ikp.kit.edu/corsika/usersguide/usersguide.pdf) 

### Create conda enviroment
```
conda create -n [ENV_NAME] --file requirements.txt
```

## Usage 
```
python reader.py --help
```

## Requirements:
 - Python 3.6
 - Numpy
 - Numba
 - Scipy
 - Astropy
