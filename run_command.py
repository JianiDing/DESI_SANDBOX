import os


def run_command(inputpath,filenames):
    '''
    :param inputpath: path that storing the data 
    :param filenames str or lists: the name/names of the hdf5 file 
    return: Running the CNN Finder Command that producing the predictions files in current directory 
    '''

    for items in filenames:
        os.system("dlacnn_anlayze_sl 1 1 DESI_MOCK2 "+str(inputpath)+str(items))
        
