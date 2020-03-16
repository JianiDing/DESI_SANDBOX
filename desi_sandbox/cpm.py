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






#defining the class for spectra
class dlametacnn:
    '''
    class for dla meta from CNN predictions

    '''

    def __init__(self):

        self.subdlaid = []

        self.z = []
        self.NHI = []
        self.dlaid = []
        self.dlac = []
        self.subz = []
        self.subNHI = []
        self.subdlac = []





    def dlametadata(self,data):




        '''
        Function for arranging dlameta from the predictions file of cnn
        :param data json file: json file of cnn predictions

        :return: class of meta data from cnn predictions, quasar redshift for detected DLAs from CNN, all mockid from detected DLAs

        '''

        mid = []
        zqso = []

        for ii in range(0, len(data)):
                # if (len(data[ii]['dlas']) > 0) or (len(data[ii]['subdlas']) > 0):
            mid.append(data[ii]['id'])
            zqso.append(data[ii]['z_qso'])
            dlanhi = []
            dlaci = []
            z_dlai = []
            subdlanhi = []
            subdlaci = []
            z_subdlai = []
            middlai = []
            midsubdlai = []
            for jj in range(0, len(data[ii]['dlas'])):
                dlanhi.append(data[ii]['dlas'][jj]['column_density'])
                dlaci.append(data[ii]['dlas'][jj]['dla_confidence'])
                z_dlai.append([data[ii]['dlas'][jj]['z_dla']])
                    # if (len(data[ii]['dlas']) > 0):
                middlai.append(data[ii]['id'])
            self.NHI.append(np.ravel(dlanhi))
            self.z.append(np.ravel(z_dlai))
            self.dlac.append(np.ravel(dlaci))  #
            self.dlaid.append(np.ravel(middlai))
                # m

            for kk in range(0, len(data[ii]['subdlas'])):
                subdlanhi.append(data[ii]['subdlas'][kk]['column_density'])
                subdlaci.append(data[ii]['subdlas'][kk]['dla_confidence'])
                z_subdlai.append(data[ii]['subdlas'][kk]['z_dla'])
                # if (len(data[ii]['subdlas']) > 0):
                midsubdlai.append(data[ii]['id'])
            self.subNHI.append(np.ravel(subdlanhi))
            self.subz.append(np.ravel(z_subdlai))
            self.subdlac.append(np.ravel(subdlaci))  #
            self.subdlaid.append(np.ravel(midsubdlai))

        return zqso,mid






def dla_output(path):
    '''
       Function for creating a dla meta object from a given directory that storing one or more predictions.json file (all the json file will be read).
       :param path for the predictions file

       :return: 1. qso redshift 2. total mockid for detected DLAs 3. meta DLA object
       '''
    totalmeta = []
    totaljson = []
    data = glob.glob(path + '/*.json')
    for ii in range(0, len(data)):
        with open(data[ii]) as f:
            totaljson.append(json.load(f))

    DLAmeta = dlametacnn()
    totalmeta = DLAmeta.dlametadata(np.concatenate(totaljson))

    return totalmeta,DLAmeta



def match_rate(dlacnnmeta, dlacnnid, dlacor, con, nhi):
    '''
    Function for matching DLAs detected by the CNN code and the DESI MOCK DLA Catalog
    :param dlacnnmeta object: meta or DLAs detected by the CNN code
    :param dlacnnid : mockid for sightlights 
    :param dlacor object: meta from DESI Mock DLA Catalog
    :param con: confidence level cut for matching
    :param nhi: NHI cut for matching
    :return:1. purity percentage  2. completeness percentage 3. meta for matched DLAs
    '''
    dlazcut = dlacor.z[np.isin(dlacor.mockid, dlacnnid)]
    dlanhicut = dlacor.NHI[np.isin(dlacor.mockid, dlacnnid)]
    dlaidcut = dlacor.mockid[np.isin(dlacor.mockid, dlacnnid)]
    
    totalmiddla = np.append(np.concatenate(dlacnnmeta.dlaid)[
                                (np.concatenate(dlacnnmeta.NHI) > nhi) & (np.concatenate(dlacnnmeta.dlac) > con)],
                            np.concatenate(dlacnnmeta.subdlaid)[
                                (np.concatenate(dlacnnmeta.subNHI) > nhi) & (np.concatenate(dlacnnmeta.subdlac) > con)])
    totalnhi = np.append(np.concatenate(dlacnnmeta.NHI)[
                             (np.concatenate(dlacnnmeta.NHI) > nhi) & (np.concatenate(dlacnnmeta.dlac) > con)],
                         np.concatenate(dlacnnmeta.subNHI)[
                             (np.concatenate(dlacnnmeta.subNHI) > nhi) & (np.concatenate(dlacnnmeta.subdlac) > con)])
    totalz = np.append(
        np.concatenate(dlacnnmeta.z)[(np.concatenate(dlacnnmeta.NHI) > nhi) & (np.concatenate(dlacnnmeta.dlac) > con)],
        np.concatenate(dlacnnmeta.subz)[
            (np.concatenate(dlacnnmeta.subNHI) > nhi) & (np.concatenate(dlacnnmeta.subdlac) > con)])
    # totalmiddla= np.concatenate(dlacnnmeta.dlaid)
    # totalnhi = np.concatenate(dlacnnmeta.NHI)
    # totalz = np.concatenate(dlacnnmeta.z)
   
    
    dlacortid = []
    dlacortz = []
    dlacortnh = []
    dlafalid = []
    dlacortmock = []
    totalmiddla = np.array(totalmiddla, dtype=int)

    count = 0
    idt = []
    print(len(totalmiddla), len(np.array(dlacor.dlaid)))
    cutid = dlaidcut[dlanhicut > nhi]
    cutnhi = dlanhicut[dlanhicut > nhi]
    cutz = dlazcut[dlanhicut > nhi]
    print (len(cutid))
    # cut2 = dla
    for ii in range(0, len(totalmiddla)):
        for jj in range(0, len(cutid)):
            if (totalmiddla[ii] == cutid[jj]) & (
                    abs(totalnhi[ii] - cutnhi[jj]) < 0.56) & (abs(totalz[ii] - cutz[jj]) < 0.015):
                count = count + 1
                #dlacortid.append(np.array(dlacor.dlaid[jj]))
                #dlacortz.append(np.array(dlacor.z[jj]))
                #dlacortnh.append(np.array(dlacor.NHI[jj]))
                #dlacortmock.append(np.array(dlacor.mockid[jj]))
    print (count)
    return count / len(totalmiddla), count / len(cutid) #dlacortmock


def confusion_matrix(TP,FP,FN,TN,version):
    '''
       Function for plotting confusion matrix for the CNN performance results
       :param TP float or int: true positive for the results
       :param FP float or int: false positive for the results
       :param FN float or int: false negative for the results
       :param TN float or int: true negative since true negative cannot be defined om this case so it is NaN
       :param: version of the mock that used to compute the confusion matrix 
       :return: plots for the confusion matrix
    '''
    plt.clf()
    cm =  [[TP, FP], [FN,TN]]
    plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Blues)
    classNames = ['Positive','Negative']
    plt.title('Confusion Matrix - '+ str(version)+' Data')
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    tick_marks = np.arange(len(classNames))
    plt.xticks(tick_marks, classNames, rotation=45)
    plt.yticks(tick_marks, classNames)
    s = [['TP','FP'], ['FN', 'NaN']]
    for i in range(2):
        for j in range(2):
            if j >= 1 & i >= 1:
                plt.text(j,i, str(s[i][j]))

            else:
                plt.text(j,i, str(s[i][j])+" = "+str(cm[i][j]))
    plt.show()



test = dla_output('/home/jding/desi_data_build/v05test1')
