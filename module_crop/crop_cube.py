# Crop a box centered on (y_pix, x_pix)




def crop_cube(x_center, y_center, cubepath, spatialcrop_pix, lam_obs, dwave):

    
    sp_crop_size = spatialcrop_pix / 2 
    
    hdul = fits.open(cubepath)
    data = hdul[1].data  
    var_data  = hdul[2].data


    y_min = int(max(0, y_center - sp_crop_size))
    y_max = int(min(data.shape[1], y_center + sp_crop_size))
    x_min = int(max(0, x_center - sp_crop_size)) 
    x_max = int(min(data.shape[2], x_center + sp_crop_size)  ) 

    sp_cropped_cube = data[:, y_min:y_max, x_min:x_max]
    var_cropped_cube = var_data[:, y_min:y_max, x_min:x_max]

    hdul[1].data = sp_cropped_cube
    hdul[2].data = var_cropped_cube
    header = hdul[1].header
    

    #  WCS spectral coordinates
    crval3 = header['CRVAL3']  # Coordinate value at the reference pixel
    cdelt3 = 1.25
    crpix3 = header['CRPIX3']  # Reference pixel index (1-indexed based on FITS standard)
    
    # converting central wavelength (lam_obs in Angstroms) to a 0-indexed pixel channel
    # Formula is pixel = ((Value - CRVAL) / CDELT) + (CRPIX - 1)
    lambda_0_pixel = ((lam_obs - crval3) / cdelt3) + (crpix3 - 1)
    
    # convert physical window width (dwave) to number of pixel channels
    deltalambda_pixel = abs(dwave / cdelt3)
    
    # Calculate integer matrix slice points using the derived pixel bounds
    num_wave = sp_cropped_cube.shape[0]
    w_min = max(0, int(round(lambda_0_pixel - deltalambda_pixel / 2)))
    w_max = min(num_wave, int(round(lambda_0_pixel + deltalambda_pixel / 2)))

    
    # slice the spatial cube array spectrally
    final_cube = sp_cropped_cube[w_min:w_max, :, :]
    var_new = var_cropped_cube[w_min:w_max, :, :] 

    
    # update FITS data matrix and modify WCS header so wavelength calibration shifts correctly
    hdul[1].data = final_cube
    hdul[2].data = var_new     
   
    # shift the spectral reference pixel to match the new cropped starting frame
    header['CRPIX3'] = crpix3 - w_min
    
    # define clean output name ,  save
    output_filename = cubepath.replace(".fits", "_fully_cropped.fits")

    hdul.writeto(output_filename, overwrite=True)
    hdul.close()

    return final_cube 


#example usage
'''
import astropy.fits as fits
from astropy.cosmology import Planck18 as COSMO
x_center = 194.22
y_center = 180.8
cubepath = "/Users/charishall/CASSI_26/J0015/J0015_vac_SUBBED.fits"
lam_obs = oII_rest * (1 + z)    # observed [OII] center, ~7945 A

spatialcrop_pix = 200

pixscale = 0.2   
vel_kms = 5000
kpc_per_arcsec  = COSMO.arcsec_per_kpc_proper(z).value
dwave =  (2 * vel_kms / c_kms) * lam_obs  # Total width in Å        ...convert velocity step to wavelength step (since delta_lam / lam = delta_v/c, then delta_lam = (delta_v/c) * lam_obs )



cropped_cube = crop_cube(x_center, y_center, cubepath, spatialcrop_pix, lam_obs, dwave)

'''