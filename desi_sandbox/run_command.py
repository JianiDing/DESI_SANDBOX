import os
import glob



def running_data(path):
    '''

    :param path: input directory for running data
    :return:
    '''
    #getting data from the directory
    data = glob.glob(path)

    #running the data
    for ii in range(33, 35):
        os.system("dlacnn_anlayze_sl 1 1 DESI_MOCK2 " + str(data[ii][36:]))



    return "end"







