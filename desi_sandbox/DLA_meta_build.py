import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy import interpolate
import h5py
from astropy.io import fits
import os

from astropy.table import Table, Column



class dlameta:
    def __init__(self,z,NHI,dlaid,mockid):


        self.z = np.array(z)
        self.NHI = np.array(NHI)
        self.dlaid = np.array(dlaid)
        self.mockid = np.array(mockid)


def reading_data(n, m, path):
    '''
    function for requesting the spectra from desi mocks
    :param resol float: dleta v for rebin
    :param zlow float: lower limit for request redshift
    :param zhigh float: higher limit for request redshift
    :return: two class object: 1.meta data and unbinned spectra 2. xspectrum object
    '''

    # define the variable

    z = []

    NHI = []
    dlaid = []
    mockid = []

    # loop through the whole directory for requesting spectra
    specf = []
    mag = []
    NHI = []
    mockmiss = []
    item1 = os.listdir(str(path))

    for k in item1[n:m]:

        item = os.listdir(str(path) + str(k))

        for j in item:

            if os.listdir(str(path) + str(k) + '/' + str(j) + '/'):
                try:
                    data = fits.open(str(path) + str(k) + '/' + str(j) + '/truth-16-' + str(
                        j) + '.fits')

                    # print (k,j)
                    zcut = ((data[1].data['z'] > 2.65) & (data[1].data['z'] < 2.70))
                    # zcut2= ((data[1].data['z'] > 2.50) & ( data[1].data['z'] < 2.55))
                    # zcut = (data[1].data['z'] > 3.8)
                    mockmiss.append(data[1].data['TARGETID'][zcut])
                    # mockmiss.append(data[1].data['TARGETID'][zcut2])
                    z.append(data[3].data['Z_DLA'])

                    dlaid.append(data[3].data['DLAID'])
                    mockid.append(data[3].data['MOCKID'])

                    NHI.append(data[3].data['NHI'])
                except (FileNotFoundError, TypeError):
                    pass

    specf = dlameta(np.concatenate(z), np.concatenate(NHI), np.concatenate(dlaid), np.concatenate(mockid))

    # newwave = np.linspace(np.min(wavetotal[0]),np.max(wavetotal[0]),reso)

    return specf, mockmiss


def dla_output_fits(totalmeta, name):

    '''
    function of building output fits file of the catalog of the CNN detected DLAs
    :param totalmeta: input catalog
    :param name: name of the output file
    :return: fits file for the catalog
    '''

    t2 = Table()
    t2['mockid'] = np.concatenate(totalmeta.dlaid)
    t2['z'] = np.concatenate(totalmeta.z)
    t2['NHI'] = np.concatenate(totalmeta.NHI)
    t2['dlac'] = np.concatenate(totalmeta.dlac)
    t2.write(name, format='fits')
    # t2.write(name, format='fits')

    return t2

