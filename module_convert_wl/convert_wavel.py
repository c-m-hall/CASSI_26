
import numpy as np
import sys
from scipy.interpolate import interp1d
from astropy.io import fits


#function to convert wair to wvac

def air2vac(wvl,precision=1.0e-12,maxiter=100):
    '''
    Convert wavelength in air to wavelength in vacuum
    Using Eq. 8 from Morton 2000
    Solution is found iteratively

    Parameters
    ----------
    wvl : float or array
        Wavelength in air
    precision : float, optional
        The preciison beyond which iteration stops;
        default is 1.0e-12
    maxiter :  int, optional
        Maximum number of iterations used; default
        is 100

    Returns
    ------
    wvl : float or array
        Wavelength in vacuum
    '''



    # First guess of wvl in vacuum is wvl in air
    wvl = np.asarray(wvl,dtype='float')
    lvac = wvl.copy()
    count = 0
    while True:
        count += 1
        s = 1.0e4/lvac
        n = (1.0+8.34254e-5+2.406147e-2/(130-s**2)+
            1.5998e-4/(38.9-s**2))
        lair = lvac/n
        dl = lair-wvl
        lvac -= dl
        s = 1.0e4/lvac
        n = (1.0+8.34254e-5+2.406147e-2/(130-s**2)+
            1.5998e-4/(38.9-s**2))
        dlair = lvac/n-wvl
        if np.abs(np.max(dlair)) < precision:
            return lvac
        if count > maxiter:
            p = np.max(dlair)
            raise(PE.PyAValError("Maximum iteration reached. Current precision is "+str(p)))



def convert_wl_main(cubename, dir):
    #fieldname = sys.argv[1]
    datapath = dir + cubename
    outpath = dir + cubename.replace('.fits', '_vac.fits')

    hdul = fits.open(datapath)
    d = hdul[1]
    v = hdul[2]

    data = d.data
    dhdr = d.header
    var = v.data
    vhdr = v.header



    Npix = dhdr['NAXIS3']
    lstep = dhdr['CD3_3']
    wave0 = np.zeros(Npix)
    steps = np.arange(0,Npix)
    wave0[:] = dhdr['CRVAL3']
    wave = wave0+lstep*steps
    wavev = air2vac(wave)
    newwave = np.arange(np.min(wavev),np.max(wavev),1.25)
    for j in range(data.shape[1]):
        print(j) 
        for i in range(data.shape[2]):
            #print(len(wavev))
            #print(len(data[:,j,i]))
            #print(len(newwave))
            data[:,j,i] = interp1d(wavev,data[:,j,i])(newwave)
            var[:,j,i] = interp1d(wavev,var[:,j,i])(newwave)

    dhdr['NAXIS3'] = len(newwave)
    vhdr['NAXIS3'] = len(newwave)
    dhdr['CRVAL3'] = np.min(newwave)
    vhdr['CRVAL3'] = np.min(newwave)

    hdul.writeto(outpath, overwrite=True)

    hdul.close()
    return outpath


#example usage 
'''
import numpy as np
import sys
from scipy.interpolate import interp1d
from astropy.io import fits

dir = '/Users/charishall/CASSI_26/'

main('J0015m0751_dr2_zap_wpsf1.corrEBV.fits', dir)
'''
