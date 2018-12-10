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
from specdb import defs as spd_defs
from specdb.build import utils as spbu
from specdb import defs as spb_defs
from astropy.table import Table




#defining the class for spectra
class spec:
    def __init__(self,ra,dec,id,filename,wvmin,wvmax,npix, date, mag,z, zerr):


        self.ra = np.array(ra)
        self.dec = np.array(dec)
        self.id = np.array(id)
        self.mag = np.array(mag)
        self.filename = filename
        self.wvmin = np.array(wvmin)
        self.wvmax = np.array(wvmax)
        self.npix = np.array(npix)
        self.date = date
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

    ra = []
    dec = []
    id = []
    filename = []
    wvmin = []
    wvmax = []
    npix = []
    date =[]
    #loop through the whole directory for requesting spectra
    specf = []
    mag = []
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
                ra.append(data[1].data['RA_TARGET'][zcut])
                dec.append(data[1].data['DEC_TARGET'][zcut])
                id.append(data[1].data['TARGETID'][zcut])
                filename.append(np.tile([k,j,j],(len(dataz[1].data['Z'][zcut]),1)))


                mag.append(data[1].data['MAG'][zcut])
                date.append(data[1].data['NIGHT'][zcut])

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


    sp = XSpectrum1D(np.vstack(np.array(wavetotal)), np.vstack(np.array(fluxt)), np.sqrt(1/np.vstack(np.array(ivar))), verbose=False)


    #ra,dec,id,filename,wvmin,wvmax,npix, date, mag,z, zerr
    binspec = ltsu.rebin_to_rest(sp, np.zeros(len(np.array(np.concatenate(z)))), resol * u.km / u.s)
    specf = spec(np.concatenate(ra), np.concatenate(dec), np.concatenate(id), np.concatenate(filename),
                 np.min(binspec.wavelength),
                 np.max(binspec.wavelength), len(binspec.wavelength), np.concatenate(date), np.concatenate(mag),
                 np.concatenate(z), np.concatenate(zerr))

    # newwave = np.linspace(np.min(wavetotal[0]),np.max(wavetotal[0]),reso)

    return specf, binspec


# starting to wrtie the hdf5 file
def hdf5_writter(reso, lowl, highl):
    '''
    :param lowl float: lower limit for request redshift
    :param highl float: higher limit for request redshift
    :param reso: delta v in km/s to rebin the data
    :return: hdf5file written
    '''
    outfil = 'desi_dla_specdb.hdf'
    hdf = h5py.File(outfil, 'w')
    #
    gdict = {}
    # creating group
    for z in np.arange(lowl, highl, 0.05):
        # hdf.create_group('z'+str(z)+'-'+str(z+0.05))
        print(z)
        spectest1 = []
        spec2 = []
        spectest1, spec2 = reading_data(reso, z, z + 0.05)
        print(np.round(z, 2), len(spectest1.z))
        # hdf_append = h5py.File('tmp'+str(z)+'.hdf', 'w')
        # hdf.create_group('z' + str(z) + '-' + str(z + 0.05))
        f = hdf.create_group('z' + str(np.round(z, 2)) + '-' + str(np.round(z + 0.05, 2)))
        npix = len(spec2.wavelength)
        data = sdb_u.init_data(npix)
        nspec = len(spectest1.z)
        print (npix)

        spec_set = hdf['z' + str(np.round(z, 2)) + '-' + str(np.round(z + 0.05, 2))].create_dataset('spec', data=data,
                                                                                                    chunks=True,
                                                                                                    maxshape=(None,),
                                                                                                    compression='gzip')
        spec_set.resize((nspec,))

        for ii in range(nspec):
            data['flux'][0][:npix] = spec2[ii].flux  # Should be flux values
            data['sig'][0][:npix] = spec2[ii].sig  # SHould be sigma values
            # print (spec[ii].sig)
            data['wave'][0][:npix] = spec2[ii].wavelength  # Should be wavelength values
            # Fill
            spec_set[ii] = data

        # hdf.copy('z'+str(z)+'-'+str(z+0.05), hdf_append['z'+str(z)+'-'+str(z+0.05)])

        # making meta data
        # hdfnew = h5py.File('z2.8_specdb_test1.hdf', 'w')
        # group = 'z2.8-2.85'
        # _ = hdfnew.create_group(group)
        group = 'z' + str(np.round(z, 2)) + '-' + str(np.round(z + 0.05, 2))
        id_key = 'DESI_ID'
        maindb, tkeys = spbu.start_maindb(id_key)

        meta = Table()
        meta['zem_GROUP'] = spectest1.z
        meta['RA_GROUP'] = spectest1.ra
        meta['DEC_GROUP'] = spectest1.dec

        meta['EPOCH'] = 2000.
        meta['sig_zem'] = 0.
        meta['flag_zem'] = np.string_('DESI')
        meta['STYPE'] = np.string_('QSO')
        # Observation
        meta['SPEC_FILE'] = np.array(spectest1.filename, dtype=float)
        meta['DATE-OBS'] = spectest1.date
        #
        meta['GROUP_ID'] = spectest1.id
        # Spectrograph
        meta['R'] = 3000.
        meta['TELESCOPE'] = np.string_('KPNO-4m')
        meta['DISPERSER'] = np.string_('ALL')
        meta['INSTR'] = np.string_('DESI')
        meta['WV_MIN'] = 3800.  # Should be the right value
        meta['WV_MAX'] = 9900.  # Should be the right value
        meta['NPIX'] = 8000  # Should be the right value
        flag_g = spbu.add_to_group_dict(group, gdict)
        maindb = spbu.add_ids(maindb, meta, flag_g, tkeys, 'DESI_ID', first=(flag_g == flag_g))
        hdf[group]['meta'] = meta
        zpri = spb_defs.z_priority()
        print (flag_g)
    spbu.write_hdf(hdf, str('test_DESI_1'), maindb, zpri, gdict, str('v0.1'),Publisher='jding')






    hdf.close()



hdf5_writter(69.08,2.8,3.0)











