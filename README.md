# sim-to-cos
Creates mock spectra for HST/COS using data from simulations

The creation of the mock spectra includes three steps:
 1. Convolution with the COS-LSF according to the required setup(grating, central wavelength and lp). __WARNING__: this version of the code does not use the wavelength dependent LSF. Rather, it uses a single lsf that is closest to the wavelength midpoint of the simulated data. 
 2. Rebinning of the simulated data to the observed spectrum.
 3. Adding poisson noise.

## Requirements

astropy~=5.1;
numpy~=1.23.3

The given versions are the tested ones. The code might work with the older versions as well.

## How to use

There are two files required as input:
#### simulated data: 
`FITS` file with columns `z`, `tau_HeII` and `sim_id` (see `examples/example_sim_data.fits`). 

The binning of the data MUST match the dispersion of the HST/COS. 
| Grating |    Dispersion     |
|---------|-------------------|
| G130M   | 9.97 mA pixel^-1  |
| G160M   | 12.23 mA pixel^-1 |
| G140L   | 80.3 mA pixel^-1  |

#### observed data: 
`FITS` file from FaintCOS with additional column `CONT` containing the fitted continuum. 

_WARNING_:`CALIB_DERED` will be used instead of `CALIB` if present.

### Arguments

The code requres six arguments:
```bash
  --sim-data PATH_TO_SIM_DATA
                        Path to the simulation data in fits format.It must have the columns "z", "tau_HeII", and "sim_id".The bin size
                        must match the COS bin size of the desired grating!
  --cos-data PATH_TO_OBS_DATA
                        Path to the observed COS spectrum reduced with FaintCOS.Additionnally, the spectrum must have the column
                        "CONT" containing fitted continuum. If data contains "CALIB_DERED" it will be used instead of "CALIB".
  --output PATH_TO_OUTPUT
                        Path to the output file. This option must end with ".fits".
  --lp LP               lifetime-poition of HST/COS
  --grating GRATING     HST/COS grating that shoud be used for mocks
  --cenwave CENWAVE     HST/COS cenwave that shoud be used for mocks
```


### Working example
```bash
python create_mock.py --sim-data=examples/example_sim_data.fits --cos-data=examples/dummy_spec.fits --lp=4 --grating=G130M --cenwave=1291 --output=examples/mocks.fits
```

