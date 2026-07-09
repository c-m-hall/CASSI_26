#PSF subtraction on MUSE datacubes containing QSO sources 


import numpy as np
from mpdaf.obj import Cube
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d






def median_rolling_window(fi, window_size=51):
    """
    Median smoothing applied on input function fi tabulated at points xi
    with a rolling window of size window_size.
    
    Parameters:
    -----------
    xi, fi     : 1d nd arrays
                 training point location, corresponding function values, and their weights 
                 Note that for numpy polyfit, weights should be set as wi=1/sigma, not as wi=1/sigma**2
                
    window_size: an odd integer, default = 5
                 must be an odd integer
                
    Returns:
    --------
    
    fi_smoothed: 1d numpy array containing smoothed values at xi
    
    """
    # first check input window_size and order with S-G requirements
    try:
        window_size = np.abs(int(window_size))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise ValueError("window_size size must be a positive odd number")
    if window_size > np.shape(fi)[0]:
        raise ValueError("window_size is too large for input function vector")
    
    whalf = (window_size-1)//2
    # first handle non-edge case where symmetric window can be used
    fip = np.copy(fi)
    for i in range(whalf, np.shape(fi)[0]-whalf):
        fid = fi[i-whalf:i+whalf+1]
        fip[i] = np.median(fid)
        
    # leaving edge pixels unsmoothed, could be handled differently if needed
    
    return fip

def interpTwoPart(wave, ratio, w):
    # smooth and interpolate through the wavelength gap with nan values


    
    mask1 = wave<5800
    mask2 = wave>6000
    mask3 = (wave>5800)&(wave<6000)
    ratio1 = savgol_filter(ratio[mask1],window_length=w,polyorder=1)
    ratio2 = savgol_filter(ratio[mask2],window_length=w,polyorder=1)
    ratio_sm = np.zeros(wave.shape)
    ratio_sm[mask1] = ratio1
    ratio_sm[mask2] = ratio2
    ratio_sm[mask3] = np.nan
    return ratio_sm

def interpLine(wv, ratio, wv1, wv2, wv3, wv4):
    # interpolate through narrow lines
    # blue window [wv1, wv2], red window [wv3, wv4]
    mask_b = (wv>wv1)&(wv<wv2)
    mask_r = (wv>wv3)&(wv<wv4)
    mask_mid = (wv>wv2)&(wv<wv3)
    med_b = np.median(ratio[mask_b])
    med_r = np.median(ratio[mask_r])
    interp = interp1d([(wv1+wv2)/2,(wv3+wv4)/2],[med_b,med_r])
    return mask_mid, interp(wv[mask_mid])

def interpWindow(line, z, dv, ddv):
    # generate interpolation window, dv in km/s
    # for lines in the line list redshifted to z

    # output blue window [wv1, wv2] ([-2dv, -dv])
    # output red window [wv3, wv4] ([+dv, +2dv])
    redline = line #*(1+z)
    dlam = dv/2.998e5*redline
    ddlam = (ddv+dv)/2.998e5*redline
    wv1, wv2 = redline-ddlam, redline-dlam
    wv3, wv4 = redline+dlam, redline+ddlam
    return wv1, wv2, wv3, wv4





#main function 

def main(cube, x, y, z, lines, template_spec, outpath):

    lines =lines # list of lines to interpolate, in rest-frame wavelength
    z = z # fiducial redshift

    # load data and calculate r from QSO center
    #keep var in case we want to use it for weighting in the future
    cube = Cube(cube)
    wave = cube.wave.coord()
    flux = cube.data.data
    var = cube.var.data
    nx, ny = flux.shape[2], flux.shape[1]
    xc, yc = x, y # center of QSO
    x, y = np.arange(0,nx), np.arange(0,ny)
    xx, yy = np.meshgrid(x, y)
    r2 = (xx-xc)**2+(yy-yc)**2
    r = np.sqrt(r2)

    # load median spectrum as QSO template
    # small amount of smoothing is applied to reduce the noise in the template

    wave_model,med_spec = np.loadtxt(template_spec,unpack=True)
    model = interp1d(np.concatenate(([4300],wave_model,[9700])), np.concatenate(([0],med_spec,[0])))(wave)

    # subtract QSO light
    rmax =25 # subtract QSO from within radius of 25 pix
    yok, xok = np.where(r<=rmax)
    for yi, xi in zip(yok, xok):
        spec = flux[:,yi,xi]
        ratio = spec/model

        ratio_sm_new = median_rolling_window(ratio, window_size=51) 
        # this is smaller than usual: due to broad [OII] and [OIII] lines from QSO with FWHM 
        # ~2000 km/s, so we want a window that won't smooth this out
        # 51 pix corresponds to ~2000 km/s, FWHM~4000 km/s
    
        flux[:,yi,xi] = spec-model*ratio_sm_new

    cube.write(outpath,  savemask='none')



#example usage 

#input list of lines to interpolate, in rest-frame wavelength
#lines = np.array([3426.863, 3727.092, 3729.875, 3869.86, 3969.591, 
                4341.692, 4687.015, 4862.683, 4960.295, 5008.240])


#main('cubepath', x_center(pixel coord), y_center(pixel_coord), z=redshift of source , lines= lines, template_spec= 'path to median template spectrum of QSO', outpath = 'output directory + filename')
