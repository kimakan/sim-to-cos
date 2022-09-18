#!/usr/bin/env python
import argparse
import os
import sys
from astropy.convolution import convolve
from astropy.table import Table, Column, vstack
from astropy.io import ascii, fits
import numpy as np

HEII_LYA = 303.7822

PATH_TO_LSF = "coslsf/"


def get_options():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sim-data', 
                        required=True,
                        dest='path_to_sim_data',
                        help='Path to the simulation data in fits format.'+
                        'It must have the columns "z", "tau_HeII", and "sim_id".'+
                        'The bin size must match the COS bin size of the desired grating!')
    parser.add_argument('--cos-data', 
                        required=True,
                        dest='path_to_obs_data',
                        help='Path to the observed COS spectrum reduced with FaintCOS.'+
                        'Additionnally, the spectrum must have the column '+
                        '"CONT" containing fitted continuum. If data contains '+
                        '"CALIB_DERED" it will be used instead of "CALIB".')
    parser.add_argument('--output', 
                        required=True,
                        dest='path_to_output',
                        help='Path to the output file. This option must end with ".fits".')
    parser.add_argument('--lp', 
                        required=True,
                        dest='lp',
                        help='lifetime-poition of HST/COS')
    parser.add_argument('--grating', 
                        required=True,
                        dest='grating',
                        help='HST/COS grating that shoud be used for mocks')
    parser.add_argument('--cenwave', 
                        required=True,
                        dest='cenwave',
                        help='HST/COS cenwave that shoud be used for mocks')

    return parser.parse_args()


def read_lsf_file(grating, lp, cenwave):
    """Reads the COS lsf file according to the setup.

    Parameters should represent a valid setup for HST/COS.

    Returns an astropy table representation of a singe lsf file.
    """

    filename = f'aa_LSFTable_{grating}_{cenwave}_LP{lp}_cn.dat'
    try:
        lsf = ascii.read(os.path.join(PATH_TO_LSF, filename), format='basic')
    except FileNotFoundError as err:
        print(f"The lsf-file for the requested setup does not exist")
        print(f"GRATING: {grating}")
        print(f"LP: {lp}")
        print(f"CENWAVE: {cenwave}")
        sys.exit()
    return lsf


def apply_lsf(wave, transmission, lsf):
    """Convolves spectrum with the cos lsf.

    In this simple version the provided transmission is convolved with the 
    COS-LSF that is closest to the mid-point of the wave-array.
    Parameters
    ----------
    wave: array-like
        Wavelength in Angstroem
    transmission: array-like
        Transmission from the simulation
    lsf: astropy-table
        Wavelength dependent COS-LSF created from the .dat file provided 
        by STSCI

    Returns
    -------
    transmission: ndarray
        Convolved transmission.
    """
    wave_midpoint = np.median(wave)
    lsf_waves = np.array([int(w) for w in lsf.colnames])
    selected_wave = lsf.colnames[np.argmin(abs(lsf_waves - wave_midpoint))]
    selected_lsf = lsf[selected_wave]
    return convolve(np.array(transmission), np.array(selected_lsf), boundary='extend')


def resample_sim_to_obs(obs_spec, sim_wave, sim_trans):
    """Resample the simulated transmission to the observed spectrum.

    Parameters
    ----------
    obs_spec: astropy table
        Observed cos spectrum reduced with FaintCOS
    sim_wave: ndarray
        Wavelength grid array from the synthetic spectrum
    sim_trans: ndarray
        Transmission from the synthetic spectrum

    Returns
    -------
    resampled_trans: ndarray
        Transmission of the synthetic spectrum sampled to the wavelength grid 
        of the observed spectrum
    """
    cos_wave = obs_spec['WAVELENGTH']
    resampled_trans = np.zeros(shape=len(cos_wave), dtype=np.float32)
    half_bin = np.mean(cos_wave[1:]-cos_wave[:-1])/2.
    for i in range(len(resampled_trans)):
        sel_trans = np.where((sim_wave>=cos_wave[i]-half_bin) & 
                             (sim_wave<cos_wave[i]+half_bin))
        if sel_trans:
            resampled_trans[i] = np.mean(sim_trans[sel_trans])
        else:
            resampled_trans[i] = 0.0

    return resampled_trans


def make_mock(obs_spec, sim_data, grating, lp, cenwave):
    mocks = []
    lsf = read_lsf_file(grating=grating, lp=lp, cenwave=cenwave)
    sim_data.add_column(Column(data=HEII_LYA*(sim_data['z']+1.), name='WAVE'))
    sim_data.add_column(Column(data=np.exp(-sim_data['tau_HeII']), name='TRANS'))
    obs_spec = obs_spec[(obs_spec['WAVELENGTH']>=min(sim_data['WAVE'])) & 
                        (obs_spec['WAVELENGTH']<=max(sim_data['WAVE']))]
    print("Creating mocks...")
    for sim_id in np.unique(sim_data['sim_id']):
        sim_data_chunk = sim_data[sim_data['sim_id']==sim_id]
        sim_data_chunk.sort('z')

        trans = apply_lsf(wave=sim_data_chunk['WAVE'],
                          transmission=sim_data_chunk['TRANS'],
                          lsf=lsf)

        trans_res = resample_sim_to_obs(obs_spec=obs_spec,
                                        sim_wave=sim_data_chunk['WAVE'],
                                        sim_trans=trans)

        bkg_err = np.random.normal(0, np.sum(obs_spec['BKG_ERR_UP']))/len(obs_spec['BKG_ERR_UP'])

        if "CALIB_DERED" in obs_spec.colnames:
            calib = obs_spec['CALIB_DERED']
        else:
            calib = obs_spec['CALIB']

        gcounts_cont = obs_spec['BACKGROUND']+bkg_err \
                        +trans_res*obs_spec['EXPTIME']*calib \
                        *obs_spec['CONT']*obs_spec['FLAT_CORR']
        gcounts_cont[np.where(gcounts_cont<0.0)] = 0.0

        gcounts = np.random.poisson(gcounts_cont)
        flux = (gcounts - obs_spec['BACKGROUND']) \
                /(obs_spec['EXPTIME']*calib*obs_spec['FLAT_CORR'])
        tmp_spec = obs_spec.copy()
        tmp_spec['GCOUNTS'] = gcounts
        tmp_spec['FLUX'] = flux
        tmp_spec.add_column(Column(trans_res, name='TRANS'))
        tmp_spec.add_column(Column(data=np.full(shape=len(flux),
                                                fill_value=sim_id,
                                                dtype=np.int32), 
                                   name='sim_id'))
        mocks.append(tmp_spec)
        sys.stdout.write("{sim_id}/{tot_sim_id}\r"\
                         .format(sim_id=sim_id,
                                 tot_sim_id=max(sim_data['sim_id'])))
        sys.stdout.flush()
    print("Finished.")
    return mocks


if __name__ == "__main__":

    args = get_options()

    obs_spec = Table(fits.getdata(args.path_to_obs_data))
    sim_data = Table(fits.getdata(args.path_to_sim_data))
    mocks = make_mock(obs_spec, sim_data, 
                      grating=args.grating, 
                      lp=args.lp, 
                      cenwave=args.cenwave)

    mock_spectra = mocks[0]
    for mock in mocks[1:]:
        mock_spectra = vstack([mock_spectra, mock])

    mock_spectra.write(args.path_to_output, format='fits', overwrite=True)

