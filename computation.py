import numpy as np
import numpy.ma as ma 
import matplotlib.pyplot as plt 
import xarray as xr
import astropy.units as u 
import astropy.constants as c 
import os 
import sys 
import exo_k as xk 
import parameterReader as pR
import argparse

def load_gcm(diagspec):
    ## we load the data and appropriate fields
    ds = xr.open_dataset(diagspec, 
                         decode_times = False)
    lon = ds.longitude.values
    lat = ds.latitude.values
    aire = ds.aire.values
    OLR3D = ds.OLR3D.values  #units is W/m2/cm-1
    wn = ds.IR_Wavenumber.values #unit is cm-1
    Rp = ds.controle.values[4] #Planetary radius (units is m)
    
    ds.close()
    
    ##let's do some time averaging  + flipping spectral axis
    OLR3d = np.flip(np.mean(OLR3D,axis=0),axis=0)
    
    ##let's create latrad and lonrad
    latrad = np.zeros((lat.size,lon.size))
    lonrad = np.zeros((lat.size,lon.size))
    latrad[:,0] = np.pi/2. + lat*np.pi/180. #this is actually co-latitude, not latitude
    lonrad[0,:] = (lon+180.)*np.pi/180. #from [-180,180] to [0:360] mais en [0:2pi]
    #let's fill the other values. 
    # I know there's a fastest way (more pythonic) to do this, just can't remember it rn
    for llo in range(lon.size):
        latrad[:,llo] = latrad[:,0]
    for lla in range(lat.size):
        lonrad[lla,:] = lonrad[0,:]
        
    ##let's create the wavelength array (units: um)
    wl = 1.e4/np.flip(wn)
    res = {
        'lon':lonrad,
        'lat':latrad,
        'Fp':OLR3d,
        'wl':wl,
        'aire':aire,
        'Rp':Rp
    }
    return res

def compute_fp(dic):
    ##getting few constants, and initialising some arrays
    nphase = dic['Fp'].shape[-1] # If you want as many phase as longitude points. SHould become a tunable parameter
    nlat = dic['Fp'].shape[1]
    nlon = dic['Fp'].shape[-1]
    nwl = dic['wl'].size
    fp = np.zeros((nphase,nwl)) #the actual planetary flux
    ls = np.arange(0,nphase)/nphase *2*np.pi #stellar longitude
    scal = np.zeros((nphase,nlat,nlon)) #visibility/scaling function
    ##let's go into terrible loops
    ## to be optimized at some points
    ## (or by a better dev than me)
    for tt in range(nphase):
        for lla in range(nlat):
            for llo in range(nlon-1):
                scal[tt,lla,llo] = np.sin(dic['lat'][lla,llo])*np.cos(dic['lon'][lla,llo]+ls[tt])*1/np.pi
                if scal[tt,lla,llo] > 0: #if the phase is visible by observers
                    fp[tt,:] +=dic['aire'][lla,llo]/dic['Rp']**2 * scal[tt,lla,llo]*dic['Fp'][:,lla,llo]
        ## flux conversion to W/m2/um
        fp[tt,:] = fp[tt,:]*1e-4*(1e4/dic['wl'])**2 #could be changed to just *= 1e4/(dic['wl']**2)
    return fp 


def stellar_flux(dic,file,bin=False):
    #add the binning down option at some point
    ## for later commits
    data = np.loadtxt(file,skiprows=1) #units are cm-1 & W/m2/cm-1
    fs = data[:,1]*1e-4*data[:,0]**2 #now units are W/m2/um
    return np.flip(fs) 


def planet_to_star(dic):
    res = np.zeros(dic['Fpcurve'].shape)
    for tt in range(res.shape[0]):
        res[tt,:] = dic['Fpcurve'][tt,:] / dic['Fs']
    
    ##now we dilute the spectra !!!
    res*= (dic['Rp']/ dic['Rs'])**2
    return res 

def saving(dic,output):
    ## I'm bored so this is gonna by .npy files
    np.save(output,dic)
    print('savet at', output)
    
def input():
    mypars = argparse.ArgumentParser(description="""Compute phase curve from GCM,
                                     using the parameters from the parameter file 
                                     given in argument of the run""",prog='phycurve')
    mypars.add_argument("-i",type=str,help='input file',action='store')
    arg =mypars.parse_args()
    initfile = vars(arg)['i']
    if initfile ==  None:
        raise argparse.ArgumentError("I need a parameter file ! Give it to me with  -i name_of_file.par")

    if '.par' not in initfile:
        raise argparse.ArgumentTypeError("""the parameter file should be passed as an argument,
                                         Give my the string of the name of the file as input,
                                         with:  -i name_of_file.par""")
    return pR.ParameterReader(initfile).dic
def main():
    ##get input parameters from the .par file 
    params = input()
    
    #attribute the init paramters
    save = params['save']
    path = params['ipath']
    specgcm =params['ifile']
    stellarfile = params['stellar_file']
    outputfile = os.path.join(params['opath'],params['ofile'])

    ## Loading data
    print("--- Loading data ---")
    dic = load_gcm(os.path.join(path,specgcm))
    dic['Rs'] = params['Rs']*u.Rsun.to(u.m) 

    #compute planetary phase curve
    print("--- Computing Fp ---")
    dic['Fpcurve'] = compute_fp(dic)
    
    ## compute stellar flux
    print("--- Computing Fs ---")
    dic['Fs'] = stellar_flux(dic,stellarfile)
    
    ##compute Fp/Fs and spectra dilution
    print("--- Computing Fp/Fs ---")
    dic['Fratio'] = planet_to_star(dic)
    
    ##Saving outputs
    if save:
        print("--- Saving outputs ---")
        saving(dic,outputfile)
        
if __name__ == '__main__':
    main()
