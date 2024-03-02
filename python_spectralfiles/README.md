# Python SpectralFile tools

This subdirectory is not part of the original SOCRATES repository. It is distributed under the 3-Clause BSD License.

### Overview

Computing k-coefficients from a line-list is very expensive, because there can be a large number of transitions for a single molecule. To speed this up, it is useful to split the process into two steps:
1. Calculating cross-sections from a line-list (at line resolution).    
2. Calculating k-coefficients from the cross-sections (at a lower resolution).    

Step 1 is by far the slowest, so it is helpful to use pre-computed opacities. These tools aim to generalise the creation of SpectralFiles from pre-computed opacity tables by converting them into a common netCDF format. These netCDF files can then be read by SOCRATES when generating SpectralFiles with the required properties.

You will need to download the source cross-sections manually. For a given database (`db`) and absorber (`ab`), the files should be placed in the directory `data/db/ab/`. 

Output files will be written to `output/`. You will need to create this directory yourself, or create a symbolic link called `output` which points to another extant location.

### Content

| Tool    | Description |
|-----------------------|-------------|
| `dace2netcdf.py`      | Convert DACE binary files into the netCDF format   | 
| `dace2xsc.py`         | Convert DACE binary files into the HITRAN xsc format   | 
| `exomol2netcdf.py`    | Convert ExoMol sigma files into the netCDF format  |
| `hitran2netcdf.py`    | Convert HITRAN xsc files into the netCDF format    |
| `wizard.py`           | Interactive wizard for generating spectral files |



### Requirements

* SOCRATES
   - Must have been compiled
   - Executables have been added to your `PATH` using `set_rad_env`
* Python (version >=3.11) and its libraries
   - numpy
   - matplotlib
   - netcdf4
   - glob
   - shutil
   - chemicals

### Online sources
* ExoMol (https://exomol.com/data/data-types/xsec/)
* DACE (https://dace.unige.ch/opacityDatabase/)
* HITRAN (https://hitran.org/xsc/)

