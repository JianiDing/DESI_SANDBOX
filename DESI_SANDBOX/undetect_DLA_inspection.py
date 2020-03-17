import matplotlib.pyplot as plt
plt.style.use(['seaborn-darkgrid'])
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy import interpolate
import h5py
from astropy.io import fits
import os
#from linetools.spectra.xspectrum import XSpectrum1D
from linetools.spectra import utils as ltsu
from astropy import units as u
from specdb.build import utils as sdb_u
from specdb import defs as spd_defs
from specdb.build import utils as spbu
from specdb import defs as spb_defs
from astropy.table import Table
import json
from linetools.spectra.xspectrum1d import XSpectrum1D
from pprint import pprint
import os
import glob
from astropy.table import Table, Column







class SPEC:
    '''
    class object for mock spectra and meta
    '''

    def __init__(self):
        self.ra = []
        self.dec = []
        self.id = []

        self.wvmin = []
        self.wvmax = []
        self.npix = []
        self.z = []
        self.zerr = []

    def spec_meta(self, ra, dec, mockid, wvmin, wvmax, npix, z, zerr ):
        self.ra = np.array(ra)
        self.dec = np.array(dec)
        self.id = np.array(mockid)

        self.wvmin = np.array(wvmin)
        self.wvmax = np.array(wvmax)
        self.npix = np.array(npix)
        self.z = np.array(z)
        self.zerr = np.array(zerr)

    def spec_meta_short(self,ra,dec,mockid,z, zerr):


        self.ra = np.array(ra)
        self.dec = np.array(dec)
        self.id = np.array(mockid)

        self.z = np.array(z)
        self.zerr = np.array(zerr)

    def __getitem__(self, sliced):
        return self.ra[sliced], self.dec[sliced], self.id[sliced], self.z[sliced], self.zerr[sliced]


def miss_id_check(path, idt, resol):
    '''
    Function for requesting spectra and meta for undetected DLAs
    :param path object: path for the spectra
    :param idt float: undetected id
    :param resol: resolution for spectr
    :return: two class object: 1. meta for spectra 2. spectra for undetected DLAs
    '''
    # define the variable
    dataz = []
    z = []
    zerr = []
    fluxt = []
    wavetotal = []
    ivar = []
    mockid = []
    ra = []
    dec = []
    mockid = []
    filename = []
    wvmin = []
    wvmax = []
    npix = []
    date = []
    # loop through the whole directory for requesting spectra
    specf = []
    mag = []
    item1 = os.listdir(str(path))
    for k in item1:
        item = os.listdir(str(path) + str(k))
        # print (item)
        for j in item:

            if os.listdir(str(path) + str(k) + '/' + str(j) + '/'):
                dataz = fits.open(
                    str(path) + str(k) + '/' + str(j) + '/truth-16-' + str(
                        j) + '.fits')
                data = fits.open(str(path) + str(k) + '/' + str(
                    j) + '/spectra-16-' + str(j) + '.fits')  # k j j
                indexcut = np.isin(data[1].data['TARGETID'], idt)
                wavecut = (np.array(data[2].data) < 5730)
                wavecut2 = ((np.array(data[7].data) > 5730) & (np.array(data[7].data) < 7560))
                wavecut3 = (np.array(data[12].data) > 7560)

                z.append(dataz[1].data['Z'][indexcut])
                zerr.append(dataz[1].data['TRUEZ'][indexcut])
                ra.append(data[1].data['TARGET_RA'][indexcut])
                dec.append(data[1].data['TARGET_DEC'][indexcut])
                mockid.append(data[1].data['TARGETID'][indexcut])

                # combining the spectra from three channels
                wavetotal.append(np.repeat(
                    [np.concatenate([data[2].data[wavecut], data[7].data[wavecut2], data[12].data[wavecut3]])],
                    len(dataz[1].data['Z']), axis=0)[indexcut])
                fluxt.append(
                    np.concatenate((data[3].data[:, wavecut], data[8].data[:, wavecut2], data[13].data[:, wavecut3]),
                                   axis=1)[indexcut])
                ivar.append(
                    np.concatenate((data[4].data[:, wavecut], data[9].data[:, wavecut2], data[14].data[:, wavecut3]),
                                   axis=1)[indexcut])

    sp = XSpectrum1D(np.vstack(np.array(wavetotal)), np.vstack(np.array(fluxt)), np.sqrt(1 / np.vstack(np.array(ivar))),
                     verbose=False)

    # ra,dec,id,filename,wvmin,wvmax,npix, date, mag,z, zerr
    binspec = ltsu.rebin_to_rest(sp, np.zeros(len(np.array(np.concatenate(z)))), resol * u.km / u.s)

    specf = SPEC.spec_meta(np.concatenate(ra), np.concatenate(dec), np.concatenate(mockid),
                 np.min(binspec.wavelength),
                 np.max(binspec.wavelength), len(binspec.wavelength),
                 np.concatenate(z), np.concatenate(zerr))

    # newwave = np.linspace(np.min(wavetotal[0]),np.max(wavetotal[0]),reso)

    return specf, binspec


def plotting(spectest, meta, n1, n2, dlacormeta, minwave, maxwave):
    '''
    plotting function for examining undetected DLA.

    :param spectest object: spectra for plotting
    :param meta float: meta for spectra
    :param n1: cloumn for plotting
    :param n2: number
    :param dlacormeta: dla meta for dla to matched
    :param minwave: maximum wavelength in plotting
    :param maxwave: minimum wavelength in plotting
    :return: two class object: plotting for undetected DLAs




    '''

    fig, ax = plt.subplots(n1, 2, figsize=(20, 18), gridspec_kw={'wspace': 0.08, 'hspace': 0.25}, sharex=False,
                           sharey=False)

    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.9, top=0.9)
    k = n2

    degrees = [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (0, 3)]
    for index in range(ax.shape[0]):
        for jndex in range(ax.shape[1]):
            dlaz = []
            nhi = []
            for ii in range(0, len(dlacormeta.mockid)):

                if dlacormeta.mockid[ii] == meta.id[k]:
                    dlaz.append(dlacormeta.z[ii])
                    nhi.append(dlacormeta.NHI[ii])
            print(dlaz)
            for ii in range(0, len(dlaz)):
                ax[index][jndex].axvline(x=(1 + dlaz[ii]) * 1215.67 / (1 + meta.z[k]), c='k',
                                         label='zabs with Nhi =' + str(nhi[ii]))
            wavecut = (np.array(spectest[k].wavelength / (1 + meta.z[k])) > 1180) & (
                        np.array(spectest[k].wavelength / (1 + meta.z[k])) < 1200)
            snr = np.mean(spectest[k].flux[wavecut] / spectest[k].sig[wavecut])

            ax[index][jndex].plot(spectest[k].wavelength / (1 + meta.z[k]), spectest[k].flux,
                                  label=str(meta.id[k]) + '/ ' + str(meta.z[k]) + '/ ' + str(snr))

            ax[index][jndex].set_xlim(minwave, maxwave)
            ax[index][jndex].set_ylim(-2, 10)
            # ax[index][jndex].axvline(rest_wave[i]*3.309,c='k',linestyle = '--',label = str(names[i]))
            # ax[index][jndex].axvline(rest_wave[i]*3.309,c='k',linestyle = '--')
            # ax[index][jndex].annotate('0.25 on axes', (0.25,4530), textcoords='data', size=20)
            ax[index][jndex].legend()
            ax[index][jndex].set_xlabel('rest frame wavelength')
            ax[index][jndex].set_ylabel('flux')
            # ax[index][jndex].set_xticklabels('{:.2f}'.format(np.array(data_indi_wave[i]/1215.67*(1+data_z[i]))))

            k = k + 1
    plt.xlim(912, 1210)
    plt.savefig('balcheck.png')

    plt.show(fig)


def snr_cal(SPEC, meta, wavemin, wavemax, cut):
    '''
    function for snr calculation

    :param SPEC object: spectra for snr calculation
    :param meta float: meta for spectra

    :param minwave: maximum wavelength in snr calculation
    :param maxwave: minimum wavelength in snr calculation
    :param cut: snr cut
    :return: snrt: 1. snr for input spectra 2. snrcutspec spectra after snr cut 3. meta for output spectra

    '''

    snrt = []
    for k in range(0, SPEC.nspec):
        wavecut = (np.array(SPEC[k].wavelength / (1 + meta.z[k])) > wavemin) & (
                    np.array(SPEC[k].wavelength / (1 + meta.z[k])) < wavemax)
        snrt.append(np.mean(SPEC[k].flux[wavecut] / SPEC[k].sig[wavecut]))
    snrcutspec = SPEC[np.where(np.array(snrt) > cut)]
    print(np.where(np.array(snrt) > cut))
    finalmeta = meta.__getitem__(np.where(np.array(snrt) > cut))

    return snrt, snrcutspec, finalmeta
