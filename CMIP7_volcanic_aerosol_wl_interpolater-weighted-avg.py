'''
Created on Aug 19, 2024

@author: Ben Johnson

code review: May Chim (Oct 21, 2024)

revised: March 04, 2025

Name of script: CMIP7_volcanic_aerosol_wl_interpolater-weighted-avg.py

Purpose:

This python script was written to aid the processing of volcanic (stratospheric) aerosol optical property data disseminated
as part of the Seventh Phase of the Couple Model Intercomparison Project (CMIP7). The script serves as an initial step
interpolating the optical properties in wavelength space on to the spectral bands of the UK Met Office Unified Model. 
The outputs of this program will be used to generate ancilary files for the UM EasyAerosol scheme 
(a subcomponent of the radiation scheme used to read in netcdf files containing climatologies of aerosol properties). 

Method: 

- Create arrays containing the waveband limits associated with an atmospheric model spectral bands

- Read the set of netcdf data files containing the CMIP7 stratospheric aerosol optical property climatologies 
  that are provided relevant to a specific set of wavelengths through the solar and terrestrial infrared spectra.
  So far this has been tested with versions 1.1.3 and 1.3.0

- Construct plank curves for the irradiance associated with solar and terrestrial radiation for use in weighted averaging

- Interpolate the optical properties from the source data onto a finer grid with many intervals within each waveband

- Integrate across each waveband to calculate the aerosol optical properties via weighted averaging. The weights are as follows:
    EXT is weighted by irradiance
    SSA is weighted by irradiance * extinction
    ASY is weighted by irradiance * extinction * SSA

- Output the interpolated aerosol optical properties to a new set of netcdf files and/or numpy files

- Option to Plot to screen a figure showing the global time mean aerosol optical properities as a function of wavelength
  for the original data and interpolated data, providing a cursory check for the outcome of the interpolation. 


External use:

This script has been shared with the anticipation that it may be useful for others involved in the preparing these input data
for use in the Seventh phase of the Coupled model Intercomparison Project (CMIP7) and is subject to the copyright below.

Please be aware that this code is shared on a good will basis and whilst is has been reviewed and tested we can not make any 
guarantees of its scientific accuracy or applicability. Users should inspect and test the outcome of running it and reach 
their own conclusion as to whether if it is the most appropriate interpolation method for their model, 
and has produced a scientifically correct / appropriate result.

Copyright 2025 UK Met Office

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this 
software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

'''

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pathlib import Path

##############################################################################
# Settings required for the weighted average method
#
# 1. Set up the plank functions with appropriate temperatures. 
# 255K is the black-body temperature of the Earth
# 5778 is representative of solar radiation
#
Trep_earth = 255.
Trep_solar = 5778.
nplank_wvs = 301             # A value of 301 seemed optimal, increasing further had no significant impact in tests
nsub_intervals_per_band = 30 # Here answers converged at around 30 and method became slower with a higher number.

##############################################################################
# Modellers will need to specify their waveband limits here or read in from a file as appropriate
#
##############################################################################
# For other models the waveband limits here will need to be altered or read in from a file as appropriate
#
def obtain_waveband_limits(wv_min_sw, wv_max_sw, wv_min_lw, wv_max_lw):
        
    # Combine sw and lw, which enable the interpolation across all wavelengths
    # rather than separately, which makes the interpolation method faster
    #
    wv_min_all = np.append(wv_min_sw, wv_min_lw)
    wv_max_all = np.append(wv_max_sw, wv_max_lw)

    # Define mid-point values of the radiation scheme wavebands
    #
    wv_mid_sw  = 0.5 * (wv_min_sw  + wv_max_sw)
    wv_mid_lw  = 0.5 * (wv_min_lw  + wv_max_lw)
    wv_mid_all = 0.5 * (wv_min_all + wv_max_all)

    # Define 2d-arrays giving the wavelength bounds for each band 
    #
    wv_bnds_all = np.column_stack((wv_min_all, wv_max_all))
    wv_bnds_sw  = np.column_stack((wv_min_sw,  wv_max_sw ))
    wv_bnds_lw  = np.column_stack((wv_min_lw,  wv_max_lw ))

    # The number of SW and LW bands
    nsw = len(wv_min_sw)
    nlw = len(wv_min_lw)

    # Specific to the UM. Some LW bands have excluded bands within them
    # If this doesn't feature in your model set exclusion_lw to False and excluded_band_lw to 0
    # just remove these lines and remove these from the arguments in the call to weighted_avg_interpolation
    #
    exclusion_lw     = [False, False, True, False, True, False, False, False, False]
    excluded_band_lw = [0,     0,      4,    0,     6,    0,     0,     0,     0]

    return wv_mid_sw, wv_mid_lw, wv_mid_all, wv_bnds_sw, wv_bnds_lw, wv_bnds_all, nsw, nlw, exclusion_lw, excluded_band_lw 

def interpolate_to_wavebands(ext_source, ssa_source, asy_source, 
                             wvs_source, 
                             wv_min_sw, wv_max_sw, 
                             wv_min_lw, wv_max_lw,
                             exclusion_lw, excluded_band_lw):
    
        
    # Set up the plank irradiance curve
    terrestrial_irrad_wavelength = np.geomspace(np.min(wv_min_lw), np.max(wv_max_lw), num=nplank_wvs)
    terrestrial_irrad            = plank_function(terrestrial_irrad_wavelength, Trep_earth)

    # Set up the solar irradiance curve
    Trep_solar = 5778.
    solar_irrad_wavelength = np.geomspace(np.min(wv_min_sw), np.max(wv_max_sw), num=nplank_wvs)
    solar_irrad            = plank_function(solar_irrad_wavelength, Trep_solar)
        
    ext_sw_interp, ssa_sw_interp, asy_sw_interp \
        = weighted_avg_interpolation(ext_source, ssa_source, asy_source, wvs_source, wv_min_sw, wv_max_sw, 
                                                solar_irrad, solar_irrad_wavelength)
        
    ext_lw_interp, ssa_lw_interp, asy_lw_interp \
        = weighted_avg_interpolation(ext_source, ssa_source, asy_source, wvs_source, wv_min_lw, wv_max_lw, 
                                                terrestrial_irrad, terrestrial_irrad_wavelength, exclusion_lw, excluded_band_lw)


    return ext_sw_interp, ssa_sw_interp, asy_sw_interp, ext_lw_interp, ssa_lw_interp, asy_lw_interp

def weighted_avg_interpolation(ext, ssa, asy, wvs_source_data, band_min, band_max, irrad, irrad_wavelength, exclusion=None, excluded_band=None):            
    
    nbands = len(band_min)
    nwvs_grid = nsub_intervals_per_band

    # Create the new arrays that will contain the optical properties interpolated onto the desired wavelengths
    shape  = np.shape(ext)
    ntimes = shape[0]
    nlats  = shape[1]
    nalts  = shape[2]

    ext_out = np.zeros((ntimes, nlats, nalts, nbands))
    ssa_out = np.zeros((ntimes, nlats, nalts, nbands))
    asy_out = np.zeros((ntimes, nlats, nalts, nbands))
    
    # Interpolation functions
    interp_func_ext = interp1d(wvs_source_data, ext, axis=3, kind='linear', bounds_error=False, fill_value="extrapolate")
    interp_func_ssa = interp1d(wvs_source_data, ssa, axis=3, kind='linear', bounds_error=False, fill_value="extrapolate")
    interp_func_asy = interp1d(wvs_source_data, asy, axis=3, kind='linear', bounds_error=False, fill_value="extrapolate")

    for iband in np.arange(0,nbands):

        print('Wavelength interpolation for band ', str(iband+1), ' of ', nbands)
        
        # Create a fine grid of wavelengths and irradiance within this band
        wvbnds_grid =  np.geomspace(band_min[iband], band_max[iband], num=nwvs_grid+1)
        wvs_grid    =  0.5 * (wvbnds_grid[1:] + wvbnds_grid[0:-1])
        irrad_grid  = np.interp(wvs_grid, irrad_wavelength, irrad)

        # wavelength increments
        dwv_grid = wvbnds_grid[1:] - wvbnds_grid[0:-1]

        # All weighted averages have this basic weighting
        basic_weighting = irrad_grid * dwv_grid
    
        # If required exclude wavelengths within the excluded band
        if exclusion != None:
            if exclusion[iband]:
                print('Excluding band ', excluded_band[iband], ' from within band ', str(iband+1))
                ind_this_band = np.where( ((wvs_grid >= band_min[iband])                    & 
                                           (wvs_grid <= band_min[excluded_band[iband]-1] ))  |
                                          ((wvs_grid >= band_max[excluded_band[iband]-1])   &
                                           (wvs_grid <= band_max[iband] )) )
                
                wvs_grid        = wvs_grid[ind_this_band]
                irrad_grid      = irrad_grid[ind_this_band]
                dwv_grid        = dwv_grid[ind_this_band]
                basic_weighting = basic_weighting[ind_this_band]

        # Interpolate optical properties onto the fine 
        ext_grid_this_band = interp_func_ext(wvs_grid)
        ssa_grid_this_band = interp_func_ssa(wvs_grid)
        asy_grid_this_band = interp_func_asy(wvs_grid)

        if np.isnan(basic_weighting).any() or np.any(basic_weighting == 0):
            raise Exception("Hey there are either Nans or zeros in the basic weighting")
            

        for ilat in np.arange(0,nlats):
            for ialt in np.arange(0,nalts):
                for itime in np.arange(0,ntimes):
                    
                    ext_grid = ext_grid_this_band[itime,ilat,ialt,:]
                    if np.nansum(ext_grid) > 0.:
                        ssa_grid = ssa_grid_this_band[itime,ilat,ialt,:]
                        asy_grid = asy_grid_this_band[itime,ilat,ialt,:]

                        # Extinction is weighted by irradiance (basic weighting)
                        ext_out[itime, ilat, ialt, iband] = weighted_avg(ext_grid, basic_weighting)
                    
                        # SSA is weighted by irradiance and extinction
                        weighting = basic_weighting * ext_grid
                        ssa_out[itime, ilat, ialt, iband] = weighted_avg(ssa_grid, weighting)
                    
                        # Asy is weighted by irradiance and scattering (extinction * ssa)
                        weighting = basic_weighting * ext_grid * ssa_grid
                        asy_out[itime, ilat, ialt, iband] = weighted_avg(asy_grid, weighting)

    return ext_out, ssa_out, asy_out

def weighted_avg(var, weighting):
    
    return np.sum(var * weighting) / np.sum(weighting)

def plank_function(wv, Trep):
    import math
    pi = np.pi
    h = 6.626e-34
    c = 3.0e+8
    k = 1.38e-23
    a = 2.0*h*pi*c**2
    b = h*c/(wv*k*Trep)
    plank_irrad = a/ ( (wv**5)*(math.e**b - 1.0) )

    return plank_irrad

def read_cmip7_source_data_v130(source_dir, filename_part, opt_prop_name, \
                                expected_dimension_keys, skip_uv_wvs=0,   \
                                first_month=0, last_month=1,              \
                                verbose=False):
    
    source_filepath = Path(source_dir) / "{}{}".format(opt_prop_name, filename_part)
    
    # Reading the source file
    print("\n### Opening the cmip7 datafile: \n", source_filepath)
    source_dataset = Dataset(source_filepath)
    
    # Check that the dimensions match what this script was set up to expect
    #
    dimensions = source_dataset.dimensions.keys()
    ikey = 0
    for key in dimensions:
        if key != expected_dimension_keys[ikey]:
            raise Exception('The dimensions of the source data did not match. I was expecting:' + expected_dimension_keys +' but got ', dimensions)
        ikey+=1

    # Put the CMIP data into local arrays
    #
    if verbose:
        print('Reading '+opt_prop_name+' data')
    this_opt_prop = source_dataset.variables[opt_prop_name][first_month-1 : last_month, :, :, skip_uv_wvs:]
    
    # Obtain the UM wavelengths
    #
    wvs_source = source_dataset.variables["wavelength"][skip_uv_wvs:]
    if verbose:
        print('The wavelengths of the source data are: \n', wvs_source)

    time      = source_dataset.variables['time'][first_month -1: last_month]
    
    return source_dataset, this_opt_prop, wvs_source, time

def nanmean_nonzero(x,ignore_zeros=False):
    if ignore_zeros:
        non_zero_ind = np.nonzero(x)    
        return np.nanmean(x[non_zero_ind])
    else:
        return np.nanmean(x)
    
def some_quick_plots(varname, var, var_sw0, var_lw0,  
                     wvs_source,
                     wv_rep_sw, wv_rep_lw,
                     sw_band_min, sw_band_max,
                     lw_band_min, lw_band_max,
                     add_irrad=False,
                     ylog=False,
                     ignore_zeros=False,
                     units='', verbose=False):

    if ylog:
        plot = plt.loglog
    else:
        plot = plt.semilogx
    
    # Look at the wavelength dependence    
    nwvs = len(wvs_source)
    nsw  = len(var_sw0[0,0,0,:])
    nlw  = len(var_lw0[0,0,0,:])
        
    gm_var     = np.empty(nwvs)
    gm_var_sw0 = np.empty(nsw)
    gm_var_lw0 = np.empty(nlw)
    
    for iwv in np.arange(0,nwvs):
        gm_var[iwv] = nanmean_nonzero(var[:,:,:,iwv], ignore_zeros=ignore_zeros)
    for iband in np.arange(0,nsw):
        gm_var_sw0[iband] = nanmean_nonzero(var_sw0[:,:,:,iband], ignore_zeros=ignore_zeros)
    
    for iband in np.arange(0,nlw):
        gm_var_lw0[iband] = nanmean_nonzero(var_lw0[:,:,:,iband], ignore_zeros=ignore_zeros)

    print('### Creating quick plot for '+varname)
    if verbose:    
        print('Global mean values of', varname, 'for SW bands: \n', gm_var_sw0)
        print('Corresponding to these mid-point SW wavelengths: \n ', wv_rep_sw)
        print('Global mean values of', varname, 'for LW bands: \n', gm_var_lw0)
        print('Corresponding to these mid-point LW wavelengths: \n ', wv_rep_lw)
            
    plt.figure()
    ax = plt.axes()
    plot(wvs_source, gm_var, marker='o', mec = 'white', mew = 1, label='source data')
    plot(wv_rep_sw, gm_var_sw0,  marker='o', mec = 'k', mew = 1, label='weighted-avg to SW bands')
    plot(wv_rep_lw, gm_var_lw0,  marker='o', mec = 'k', mew = 1, label='weighted-avg to LW bands')
    plt.title('Evaluating wavelength interpolation for '+varname)
    plt.ylabel(' '.join(['Global mean', varname, units]))
    plt.xlabel('Wavelength (m)')

    if add_irrad:
        plank_wvs    = np.geomspace(np.min(lw_band_min), np.max(lw_band_max), num=nplank_wvs)
        plank_irrad0 = plank_function(plank_wvs, Trep_earth)
        plank_irrad  = plank_irrad0 * 1.e-7 / np.max(plank_irrad0)
        
        solar_wvs    = np.geomspace(np.min(sw_band_min), np.max(sw_band_max), num=nplank_wvs)
        solar_irrad0 = plank_function(solar_wvs, Trep_solar)
        solar_irrad  = plank_irrad0 * 1.e-7 / np.max(solar_irrad0)

        solar_irrad = solar_irrad0 * 1.e-7 / np.max(solar_irrad0) 
        plot(solar_wvs, solar_irrad, color='k', label='Plank curve at '+str(Trep_solar)+'K (arb units)')
        plot(plank_wvs, plank_irrad, color='k', label='Plank curve at '+str(Trep_earth)+'K (arb units)')
        
    plt.legend()

# This function writes out the interpolated data to a netcdf file
# (just a single optical property per file) with the same dimensions
# and attributes as the original source file, except for the wavelength
# dimension and possible the time dimension
#           
def write_interpolated_netcdf(output_dir, interpolated_filename, opt_prop_name, \
                                   source_dataset, opt_prop_interp, wvs_source, \
                                   wv_mid, wv_bnds, time):

    # Define new filepath for file to write out
    #
    new_interpolated_ncfile = Path(output_dir) / interpolated_filename
        
    # Create a new NetCDF file
    #
    new_dataset = Dataset(new_interpolated_ncfile, 'w', format='NETCDF4')
    
    # Copy dimensions, adjusting the time dimension
    #
    for name, dimension in source_dataset.dimensions.items():
        if name == 'time':
            new_dataset.createDimension(name, len(time))
        elif name == 'wavelength':
            new_dataset.createDimension(name, len(wv_mid))
        else:
            new_dataset.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
    
    # Copy variables and their attributes
    for name, variable in source_dataset.variables.items():
        if name == 'wavelength':
            new_var = new_dataset.createVariable(name, variable.datatype, ('wavelength',))
            new_var[:] = wv_mid
            new_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            
        elif name == opt_prop_name:
            new_var = new_dataset.createVariable(name, variable.datatype, ('time', 'lat', 'height', 'wavelength'))
            new_var[:] = opt_prop_interp
            new_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            
        elif name == 'time':
            new_var = new_dataset.createVariable(name, variable.datatype, ('time',))
            new_var[:] = time
            new_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            
        elif name == 'wavelength_bnds':
            new_var = new_dataset.createVariable(name, variable.datatype, ('wavelength','bnds'))
            new_var[:] = wv_bnds
            new_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            
        else:
            new_var = new_dataset.createVariable(name, variable.datatype, variable.dimensions)
            new_var.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            new_var[:] = variable[:]
    
    # Copy global attributes
    new_dataset.setncatts({k: source_dataset.getncattr(k) for k in source_dataset.ncattrs()})
    
    # Close the new dataset
    print('### Writing interpolated data to: \n', new_interpolated_ncfile)
    if verbose:
        print('For the mid-point wavelengths of: \n', wv_mid)
        print('Corresponding to wavelength bounds of: \n', wv_bnds)
        
    new_dataset.close()

if __name__ == '__main__':

    ###############################################################################

    # START OF USER SPECIFICATIONS

    ###############################################################################
    # 1. Specific directories and filepaths for input4MIPs monthly volcanic forcing 
    # data and outputs
    ###############################################################################
    #
    source_dir   = '/your_base_directory/v1.3.0'
    
    # Data version. This could be used in the directory structure below source_dir 
    # data_version = 'v20250122'
    
    # Filenames minus the 3-charcter string at the start corresponding to the optical property name 
    #
    filename_part = '_input4MIPs_aerosolProperties_CMIP_UOEXETER-CMIP-1-3-0_gnz_185001-202112-clim.nc'

    # Output file for interpolated data
    #
    output_dir = 'your_output_dir'
    
    interpolated_filename_part  = 'CMIP7_interpolated'

    ###############################################################################
    # 2. Specify the wavelength limits of the radiation scheme bands
    ###############################################################################
        
    # These are the wavelength boundaries of the Met Office Unified Model shortwave 
    # radiation scheme, as of Global Atmosphere configuration GA7/8/9.
    #
    # Units are meters.
    
    wv_min_sw = np.array([2.0e-07, 3.2e-07, 5.05e-07, 6.9e-07, 1.19e-06, 2.38e-06])
    wv_max_sw = np.array([3.2e-07, 5.05e-07, 6.9e-07, 1.19e-06, 2.38e-06, 1.0e-05])
        
    # These are the wavelength boundaries of the Met Office Unified Model shortwave 
    # radiation scheme, as of Global Atmosphere configuration GA7/8/9
    #
    # Units are meters.
    
    wv_min_lw = np.array([25.0e-06, 18.182e-06, 12.500e-06, 13.333e-06, 8.333e-06, 8.929e-06, 7.519e-06, 6.667e-06, 3.339e-06]) 
    wv_max_lw = np.array([100.0e-06,  25.0e-06, 18.182e-06, 16.949e-06, 12.50e-06, 1.010e-05, 8.333e-06, 7.519e-06, 6.667e-06]) 

    # Users may want to avoid reading data from the source file corresponding to the first wavelength
    # if that falls outside the spectral range of the radiation scheme. 
    # If so set skip_uv_wvs = 1, other wise set to 0.
    #
    skip_uv_wvs = 0
    
    ###############################################################################
    # 3. Specify which months of the time series do you want to read 
    ###############################################################################
    #
    # As of version 1.1.3 and 1.3.0 there were 3288 months in the file 
    # covering the time period 1750 - 2023 (274 years * 12 = 3288).
    #
    # Therefore users may want to use the following setting for first_month and last_month:
    #
    # ****** !!! These will need changing if future releases of the input files 
    #******* !!! cover a different date range
    #
    Jan_1750 = 1     # First month of timeseries for any simulation starting in 1750
    Jan_1850 = 1201  # For historical simulations starting in 1850 
    Dec_2021 = 3264  # To include data through to the end of the 2021  
    Dec_2023 = 3288  # To include data through to the end of the timeseries
                     # or could set last_month to -1 to indicate reading through
                     # to final month of the timeseries   
    #
    # If reading from the 12-month historical average climatology set to 1 and 12
    first_month = 1
    last_month  = 12

    ###############################################################################
    # 4. Specify what you would like this script to do
    ###############################################################################
    #
    # Would like the interpolated data to be written out
    # to either a netcdf or numpy file (or both)? 
    #
    save_to_netcdf_file = True
    save_to_numpy_file  = False
    
    # Would you also like quick plots of the interpolated data plotted to screen?
    #
    show_quick_plots    = True
    
    # Would you  like to use the verbose option with a higher level of standard output?
    #
    verbose             = True

    ###############################################################################

    # END OF USER SPECIFICATIONS

    ###############################################################################
    

    ###############################################################################
    # This obtains the wavelength limits of the SW and LW bands of the model that 
    # we need to interpolate to
    #
    wv_mid_sw, wv_mid_lw, \
    wv_mid_all, wv_bnds_sw, \
    wv_bnds_lw, wv_bnds_all, \
    nsw, nlw, \
    exclusion_lw, excluded_band_lw \
        = obtain_waveband_limits(wv_min_sw, wv_max_sw, wv_min_lw, wv_max_lw)
    
    
    ###############################################################################
    # This lists the optical properties names for the loop below and creating the full 
    # filenames of source netcdfs
    #
    opt_prop_names = ['ext', 'ssa', 'asy']

    ###############################################################################
    # As of version 1.1.3 the optical properties have been given with the dimensions listed below.
    # The interpolation method will not work and is set to crash if the order of dimensions 
    # does not match this:
    #
    expected_dimension_keys = ['time', 'lat', 'height', 'wavelength', 'bnds', 'nv']

    ###############################################################################
    # Information only for quick check plotting function below
    #
    ylog           = [True, False, False]
    units          = ['(m$^{-1}$)', '', '']
    iopt           = 0

    ###############################################################################    
    # 
    # Get the optical properties from the source files
    #   
    source_dataset_ext, ext_source, wvs_source, time = \
                read_cmip7_source_data_v130(source_dir, filename_part, \
                                            'ext', expected_dimension_keys, \
                                            skip_uv_wvs=skip_uv_wvs, \
                                            first_month=first_month, last_month=last_month, \
                                            verbose=verbose)

    source_dataset_ssa, ssa_source, wvs_source, time = \
                read_cmip7_source_data_v130(source_dir, filename_part, \
                                            'ssa', expected_dimension_keys, \
                                            skip_uv_wvs=skip_uv_wvs, \
                                            first_month=first_month, last_month=last_month, \
                                            verbose=verbose)

    source_dataset_asy, asy_source, wvs_source, time = \
                read_cmip7_source_data_v130(source_dir, filename_part, \
                                            'asy', expected_dimension_keys, \
                                            skip_uv_wvs=skip_uv_wvs, \
                                            first_month=first_month, last_month=last_month, \
                                            verbose=verbose)
               
    # Interpolate on to the model wavebands using chosen method
    ext_sw_interp, \
    ssa_sw_interp, \
    asy_sw_interp, \
    ext_lw_interp, \
    ssa_lw_interp, \
    asy_lw_interp, \
        = interpolate_to_wavebands(ext_source, ssa_source, asy_source, wvs_source, 
                                   wv_min_sw, wv_max_sw, wv_min_lw, wv_max_lw, 
                                   exclusion_lw, excluded_band_lw)
                             
    opt_prop_sw_interp_list = [ext_sw_interp, ssa_sw_interp, asy_sw_interp]
    opt_prop_lw_interp_list = [ext_lw_interp, ssa_lw_interp, asy_lw_interp]
    source_datasets         = [source_dataset_ext, source_dataset_ssa, source_dataset_asy]
    
    # Loop over the three optical properties
    #
    for (opt_prop_name, opt_prop_sw_interp, opt_prop_lw_interp, source_dataset) in zip(opt_prop_names, opt_prop_sw_interp_list, opt_prop_lw_interp_list, source_datasets):
        
        if save_to_netcdf_file:
            
            # Write the interpolated data out to a new netcdf file that will have
            # the same attributes and dimensions of the source file, except for the new
            # set of wavelengths and possibly a subset of the time
            #
 
            # Writing out for the SW
            #
            sw_interpolated_filename = '{}{}{}{}'.format(opt_prop_name, "_sw_", interpolated_filename_part, ".nc")
            write_interpolated_netcdf( output_dir, sw_interpolated_filename, opt_prop_name, \
                                       source_dataset, opt_prop_sw_interp, wvs_source,  \
                                       wv_mid_sw, wv_bnds_sw, time)
    
            # Write out for the LW
            #
            lw_interpolated_filename = '{}{}{}{}'.format(opt_prop_name, "_lw_", interpolated_filename_part, ".nc")
            write_interpolated_netcdf( output_dir, lw_interpolated_filename, opt_prop_name, \
                                       source_dataset, opt_prop_lw_interp, wvs_source,  \
                                       wv_mid_lw, wv_bnds_lw, time)
            
        if save_to_numpy_file:
            
            # Save the interpolated data also to a numpy file if desired
            #
            output_filepath_sw = Path(output_dir) / '{}{}{}{}'.format(opt_prop_name, "_sw_", interpolated_filename_part, ".npz")
            output_filepath_lw = Path(output_dir) / '{}{}{}{}'.format(opt_prop_name, "_lw_", interpolated_filename_part, ".npz")

            print('### Writing interpolated data to numpy arrays in: \n', output_filepath_sw, '\n', output_filepath_lw)

            np.savez(output_filepath_sw, \
                     opt_prop_sw_interp=opt_prop_sw_interp)

            np.savez(output_filepath_lw, \
                     opt_prop_sw_interp=opt_prop_lw_interp)

        # close the source dataset
        # 
        source_dataset.close()

    ###############################################################################
    # If requested quick this will visualise the original and interpolated data
    # for the selected field
    #
    if show_quick_plots:
        #
        some_quick_plots('Extinction', ext_source, ext_sw_interp, ext_lw_interp, 
                         wvs_source,
                         wv_mid_sw, wv_mid_lw,
                         wv_min_sw, wv_max_sw, 
                         wv_min_lw, wv_max_lw, 
                         add_irrad=True,
                         ylog=True,
                         ignore_zeros=True,
                         units='/km')
        
        some_quick_plots('SSA', ssa_source, ssa_sw_interp, ssa_lw_interp, 
                         wvs_source,
                         wv_mid_sw, wv_mid_lw,
                         wv_min_sw, wv_max_sw, 
                         wv_min_lw, wv_max_lw, 
                         ignore_zeros=True)
            
        some_quick_plots('ASY', asy_source, asy_sw_interp, asy_lw_interp, 
                         wvs_source,
                         wv_mid_sw, wv_mid_lw,
                         wv_min_sw, wv_max_sw, 
                         wv_min_lw, wv_max_lw, 
                         ignore_zeros=True)
        plt.show()

	###############################################################################
    # END OF PROGRAM
    #        
    print('All done')
