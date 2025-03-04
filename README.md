# CMIP7_volcanic_aerosol_forcing: Wavelength interpolation tools
This repository provides two wavelength interpolation tools for processing stratospheric aerosol optical property data used to prescribe volcanic forcing during CMIP7 and may be of benefit to other climate modellers participating in CMIP7. 

CMIP7_volcanic_aerosol_wl_interpolator-midpoint.py

CMIP7_volcanic_aerosol_wl_interpolator-weighted-avg.py

Both of these offer a method for interpolating the CMIP7 volcanic aerosol optical properties across wavelength space and output either 
another netcdf or numpy file with the data on the wavelengths intervals corresponding to a model's radiation scheme. Users may wish to consider the differences between the two methods and if either will be of assistance it is strongly recommended to check the output of the interpolation to assess the suitability for their model / application. 

The midpoint version is simple and fast to run. It performs a straight-forward linear interpolation to the midpoint
wavelengths that are specified, corresponding to a radiation scheme spectral bands. The principal downside of this method is that it can lead
to unrepresentative optical properties when applied to for wide wavebands that cover sections of the solar or terrestrial
spectrum with strong gradients in optical properties and/or irradiance with wavelength.

The weighted-avg version is a more complex version that will require greater CPU and memory but offers
potentially more accurate results, particularly when applied to wide spectral wavebands where the optical
properties and solar or terrestrial radiation have strong gradients in wavelength space. This method calculates weighted averages
of the optical properties across the wavelength limits of a spectral band by splitting the band into smaller intervals and itegrating with appropriate weightings:

  Extinction:                   Weighted by irradiance
  
  Single-scattering albedo:     Weighted by irradiance * extincton
  
  Asymmetry factor:             Weighted by irradiance * extinction * single-scattering albeod

Please be aware that these tools are shared on a good will basis and whilst they have been reviewed and tested we can not make any 
guarantees of the scientific accuracy or applicability. As noted above, users should inspect and test the outcome of running it and reaching 
their own conclusion as to whether if it is the most appropriate interpolation method for their model, and has produced a scientifically correct / appropriate result. 

If either of these tools seem potentially useful please take a working copy and modify the main section of the script to specify the 
wavelengths of your model and input / output directories etc. 

See comments within code for further details. 
