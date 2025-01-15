# CMIP7_volcanic_aerosol_forcing: Wavelength interpolation tool
This repository provides a tool for processing stratospheric aerosol optical property data used to prescribe volcanic forcing during CMIP7 and may be of benefit to other climate modellers participating in CMIP7. 

CMIP7_volcanic_aerosol_wl_interpolator-midpoint.py

This interpolates the CMIP7 volcanic aerosol optical properties across wavelength space and outputs either 
another netcdf or numpy file with the data on the wavelengths intervals corresponding to your model's radiation scheme. 

If this seems potentially useful please take a working copy and modify the main section of the script to specify the 
wavelengths of your model and input / output directories etc. 

See comments within the code for details. 
