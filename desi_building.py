import matplotlib.pyplot as plt
plt.style.use(['seaborn-darkgrid'])
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy import interpolate
import h5py
from astropy.io import fits
import os
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import utils as ltsu
from astropy import units as u
from specdb.build import utils as sdb_u



#defining the class for spectra
class spec:
    def __init__(self, wave, flux, ivar, z, zerr):
        self.wave = np.array(wave)
        self.flux = np.array(flux)
        self.ivar = np.array(ivar)


        self.z = np.array(z)
        self.zerr = np.array(zerr)


def reading_data(resol, zlow, zhigh):
    '''
    function for requesting the spectra from desi mocks
    :param resol float: dleta v for rebin
    :param zlow float: lower limit for request redshift
    :param zhigh float: higher limit for request redshift
    :return: two class object: 1.meta data and unbinned spectra 2. xspectrum object
    '''

    #define the variable
    dataz = []
    z = []
    zerr = []
    fluxt = []
    wavetotal = []
    ivar = []
    mask = []

    #loop through the whole directory for requesting spectra
    specf = []
    item1 = os.listdir("/data/desi_mock/v03/quick-0.1/quick-0.1/spectra-16/")
    for k in item1:
        item = os.listdir("/data/desi_mock/v03/quick-0.1/quick-0.1/spectra-16/" + str(k))
        # print (item)
        for j in item:

            if os.listdir('/data/desi_mock/v03/quick-0.1/quick-0.1/spectra-16/' + str(k) + '/' + str(j) + '/'):
                dataz = fits.open(
                    '/data/desi_mock/v03/quick-0.1/quick-0.1/spectra-16/' + str(k) + '/' + str(j) + '/zbest-16-' + str(
                        j) + '.fits')
                data = fits.open('/data/desi_mock/v03/quick-0.1/quick-0.1/spectra-16/' + str(k) + '/' + str(
                    j) + '/spectra-16-' + str(j) + '.fits')  # k j j
                zcut = (zlow < np.array(dataz[1].data['Z'])) & (np.array(dataz[1].data['Z']) < zhigh)
                wavecut = (np.array(data[2].data) < 5730)
                wavecut2 = ((np.array(data[7].data) > 5730) & (np.array(data[7].data) < 7560))
                wavecut3 = (np.array(data[12].data) > 7560)
                z.append(dataz[1].data['Z'][zcut])
                zerr.append(dataz[1].data['ZERR'][zcut])

                #combining the spectra from three channels
                wavetotal.append(np.repeat(
                    [np.concatenate([data[2].data[wavecut], data[7].data[wavecut2], data[12].data[wavecut3]])],
                    len(dataz[1].data['Z']), axis=0)[zcut])
                fluxt.append(
                    np.concatenate((data[3].data[:, wavecut], data[8].data[:, wavecut2], data[13].data[:, wavecut3]),
                                   axis=1)[zcut])
                ivar.append(
                    np.concatenate((data[4].data[:, wavecut], data[9].data[:, wavecut2], data[14].data[:, wavecut3]),
                                   axis=1)[zcut])


    specf = spec(np.vstack(np.array(wavetotal)), np.vstack(np.array(fluxt)), np.vstack(np.array(ivar)), np.concatenate(z), np.concatenate(zerr))
    sp = XSpectrum1D(specf.wave, specf.flux, specf.ivar, verbose=False)

    binspec = ltsu.rebin_to_rest(sp, np.zeros(len(np.array(specf.wave))), resol * u.km / u.s)
    # newwave = np.linspace(np.min(wavetotal[0]),np.max(wavetotal[0]),reso)

    return specf, binspec





#print (len(spectest1.wave))

# starting to wrtie the hdf5 file
def hdf5_writter(reso,lowl,highl):
    '''
    :param lowl float: lower limit for request redshift
    :param highl float: higher limit for request redshift
    :param reso: delta v in km/s to rebin the data
    :return: hdf5file written
    '''
    outfil = 'desi_dla_1.hdf'
    hdf = h5py.File(outfil, 'w')
    #

    # creating group
    for z in np.arange(lowl,highl,0.05):
        #hdf.create_group('z'+str(z)+'-'+str(z+0.05))
        print (z)
        spectest1 = []
        spec2 = []
        spectest1, spec2 = reading_data(reso, z, z+0.05)
        print (np.round(z,2),len(spectest1.wave))
        #hdf_append = h5py.File('tmp'+str(z)+'.hdf', 'w')
        #hdf.create_group('z' + str(z) + '-' + str(z + 0.05))
        f = hdf.create_group('z' + str(np.round(z,2)) + '-' + str(np.round(z+0.05,2)))
        npix = len(spec2.wavelength)
        data = sdb_u.init_data(npix)
        nspec = len(spectest1.wave)

        spec_set = hdf['z' + str(np.round(z,2)) + '-' + str(np.round(z+0.05,2))].create_dataset('spec', data=data, chunks=True,maxshape=(None,), compression='gzip')
        spec_set.resize((nspec,))

        spec_setz = hdf['z' + str(np.round(z, 2)) + '-' + str(np.round(z + 0.05, 2))].create_dataset('zem', data=spectest1.z,chunks=True,maxshape=(None,),compression='gzip')
        spec_setz.resize((nspec,))
        #creacting z attribute for dataset but cannot exceed 20000 columns
        #spec_set.attrs['zem'] = spectest1.z
        for ii in range(nspec):




            data['flux'][0][:npix] = spec2[ii].flux  # Should be flux values
            data['sig'][0][:npix] = spec2[ii].sig  # SHould be sigma values
            #print (spec[ii].sig)
            data['wave'][0][:npix] = spec2[ii].wavelength  # Should be wavelength values
        # Fill
            spec_set[ii] = data




        #hdf.copy('z'+str(z)+'-'+str(z+0.05), hdf_append['z'+str(z)+'-'+str(z+0.05)])

    hdf.close()


#creating hdf5 file by using a delta v for 69.08 km/s for a redshift range from 1.9-2.0
#hdf5_writter(69.08,2.8,3.0)






#checking hdf5
#hdf_new = h5py.File('desi_dla_test.hdf', 'r')
#print (list(hdf_new.keys()))
#print (hdf_new['z2.9-2.95']['zem'].value)

#append testing

#hdf_append = h5py.File('tmp2.hdf', 'w')
#hdf_append.create_group('z1.9-2.0')
#hdf.copy('z1.9-2.0/spec', hdf_append['z1.9-2.0'])

