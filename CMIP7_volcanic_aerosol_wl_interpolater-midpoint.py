'''
Created on Aug 19, 2024

@author: Ben Johnson

code review: May Chim (Oct 21, 2024)

Name of script: CMIP7_volcanic_aerosol_wl_interpolater-midpoint.py

Purpose:

This python script was written to aid the processing of volcanic (stratospheric) aerosol optical property data disseminated
as part of the Seventh Phase of the Couple Model Intercomparison Project (CMIP7). The script serves as an initial step
interpolating the optical properties in wavelength space on to the spectral bands of the UK Met Office Unified Model. 
The outputs of this program will be used to generate ancilary files for the UM EasyAerosol scheme 
(a subcomponent of the radiation scheme used to read in netcdf files containing climatologies of aerosol properties). 

Method: 

- Create arrays containing the waveband limits and midpoint wavelengths associated with an atmospheric model spectral bands

- Read the set of netcdf data files containing the CMIP7 stratospheric aerosol optical property climatologies (as at version 1.1.3),
  that are provided relevant to a specific set of wavelengths through the solar and terrestrial infrared spectra.
  
- Interpolate the aerosol optical properties linearly from the input data onto the mid-point of the established spectral bands

- Output the interpolated aerosol optical properties to a new set of netcdf files and/or numpy files

- Option to Plot to screen a figure showing the global time mean aerosol optical properities as a function of wavelength
  for the original data and interpolated data, providing a cursory check for the outcome of the interpolation. 


External use:

This script has been shared with the anticipation that it may be useful for others involved in the preparing these input data
for use in the Seventh phase of the Coupled model Intercomparison Project (CMIP7) and is subject to the copyright below.

BSD 3-Clause License:

Copyright (c) 2025, mo-benjohnson (Met Office, UK)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from pathlib import Path

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

    return wv_mid_sw, wv_mid_lw, wv_mid_all, wv_bnds_sw, wv_bnds_lw, wv_bnds_all, nsw, nlw 

def mid_point_interpolation(opt_prop_source, wvs_source, band_mid, verbose=False):                    

    print('### Starting mid-point interpolation')

    if verbose:
        print('The full list of midpoint wavelengths including both SW and LW bands that we are interpolating to are: \n', band_mid)
        
    # Remarkably simple function thanks to scipy
    #
    interp_func = interp1d(wvs_source, opt_prop_source, axis=-1, kind='linear', bounds_error=False, fill_value="extrapolate")
    #
    opt_prop_interp = interp_func(band_mid)
    
    return opt_prop_interp

def read_cmip7_source_data_v130(source_dir, data_version, filename_part,    \
                                opt_prop_name, expected_dimension_keys,     \
                                skip_uv_wvs=0, first_month=0, last_month=1, \
                                verbose=False):
    
    # The first wavelength in the source data (0.1 microns) may be well outside the range 
    # of the SW spectral bands so this can be skipped if skip_uv_wvs is set to 1
    #
    source_filepath = Path(source_dir) / data_version / "{}{}".format(opt_prop_name, filename_part)
    
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
    time_bnds = source_dataset.variables['time_bnds'][first_month -1 : last_month]
    
    return source_dataset, this_opt_prop, wvs_source, time, time_bnds 

# This little function is only used in quick_plots
#
def nanmean_nonzero(x,ignore_zeros=False):
    if ignore_zeros:
        non_zero_ind = np.nonzero(x)    
        return np.nanmean(x[non_zero_ind])
    else:
        return np.nanmean(x)

# This function shows the global mean interpolated and original values in a plot to screen   
#
def quick_plots(varname, var, var_sw0, var_lw0,  
                     wvs_source,
                     wv_rep_sw, wv_rep_lw,
                     nsw, nlw,
                     ylog=False,
                     ignore_zeros=False,
                     units='', verbose=False):

    if ylog:
        plot = plt.loglog
    else:
        plot = plt.semilogx
    
    # Look at the wavelength dependence    
    nwvs = len(wvs_source)
   
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
    plot(wv_rep_sw, gm_var_sw0,  marker='o', mec = 'k', mew = 1, label='mid-point interp SW bands')
    plot(wv_rep_lw, gm_var_lw0,  marker='o', mec = 'k', mew = 1, label='mid-point interp LW bands')
    plt.title('Evaluating wavelength interpolation for '+varname)
    plt.ylabel(' '.join(['Global mean', varname, units]))
    plt.xlabel('Wavelength (m)')
    
    plt.legend()

# This function writes out the interpolated data to a netcdf file
# (just a single optical property per file) with the same dimensions
# and attributes as the original source file, except for the wavelength
# dimension and possible the time dimension
#                                 
def write_interpolated_netcdf(output_dir, interpolated_filename, opt_prop_name, \
                                   source_dataset, opt_prop_interp, wvs_source, \
                                   wv_mid, wv_bnds, time, time_bnds):

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
            
        elif name == 'time_bnds':
            new_var = new_dataset.createVariable(name, variable.datatype, ('time','bnds'))
            new_var[:] = time_bnds
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

    # Data version. This should be used in the directory structure below source_dir 
    data_version = 'v20250122'
    
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

    # Users may want to avoid reading data from the source file corresponding to the first wavelength (0.1 microns) 
    # if that falls well outside the spectral range of the radiation scheme. 
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
                     # (or in this case could set last_month to -1) 
    #
    # Assuming we'd like data covering the period 1850 - 2021
    first_month = Jan_1850
    last_month  = Dec_2021

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
    wv_mid_sw, wv_mid_lw, wv_mid_all,       \
    wv_bnds_sw, wv_bnds_lw, wv_bnds_all,    \
    nsw, nlw =                              \
        obtain_waveband_limits(wv_min_sw, wv_max_sw, wv_min_lw, wv_max_lw)
    
    
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
    # Loop over the three optical properties
    #
    for opt_prop_name in opt_prop_names:

        # Read the CMIP7 source data file to access the information we need
        # for aerosol optical properties interpolation
        #
        
        source_dataset, opt_prop_source, wvs_source, time, time_bnds = \
                                                                            \
                    read_cmip7_source_data_v130(source_dir, data_version, filename_part, \
                                           opt_prop_name, expected_dimension_keys, skip_uv_wvs=skip_uv_wvs, \
                                           first_month=first_month, last_month=last_month, \
                                           verbose=verbose)

        # Perform linear interpolation from the source data to the mid-point of model wavebands
        #
        opt_prop_all_interp = mid_point_interpolation(opt_prop_source, wvs_source, wv_mid_all, verbose=verbose)    

        # Separate out into SW and LW components, assuming model radiation scheme treat these
        # separately / want the three input files (ext, ssa, asy) for each spectrum
        #
        opt_prop_sw_interp = opt_prop_all_interp[:,:,:,0:nsw]
        opt_prop_lw_interp = opt_prop_all_interp[:,:,:,nsw:nsw+nlw]
        
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
                                       wv_mid_sw, wv_bnds_sw, time, time_bnds)
    
            # Write out for the LW
            #
            lw_interpolated_filename = '{}{}{}{}'.format(opt_prop_name, "_lw_", interpolated_filename_part, ".nc")
            write_interpolated_netcdf( output_dir, lw_interpolated_filename, opt_prop_name, \
                                       source_dataset, opt_prop_lw_interp, wvs_source,  \
                                       wv_mid_lw, wv_bnds_lw, time, time_bnds)
            
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
            quick_plots(opt_prop_name, opt_prop_source, opt_prop_sw_interp, opt_prop_lw_interp, 
                             wvs_source,
                             wv_mid_sw, wv_mid_lw,
                             nsw, nlw,
                             ylog=ylog[iopt],
                             ignore_zeros=True,
                             units=units[iopt], 
                             verbose=verbose)
            iopt+=1
            
    # If requested plot to screen
    #
    if show_quick_plots:
        plt.show()

    ###############################################################################
    # END OF PROGRAM
    #
    print('All done')
