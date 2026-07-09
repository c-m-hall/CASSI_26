# CASSI_26 -> stacking MUSE datacubes



This repo contains the code for preparing MUSE datacubes for stacking. 

The pipeline is as follows: 
  - download the cube
  - convert the wavelength from air to vaccuum
  - use PSF subtraction to remove continuum flux at the core of the source
  - crop the cube in the spatial and spectral directions
  - resample/interpolate onto a new grid (while conserving flux)
  - build a sky spectrum mask (using median absolute deviation), and mask residuals

