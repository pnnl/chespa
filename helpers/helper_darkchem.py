# -*- coding: utf-8 -*-
"""
Helper file to get DarkChem latent space vectors used in paper

Note: DarkChem must already be downloaded.
For more info, see https://github.com/pnnl/darkchem

Once downloaded, set DARKCHEM_PATH variable

"""

## Imports ##
import pandas as pd
import os
import numpy as np


## Globals ##
# save path to trained protonated weights, relative to the path of the code
DARKCHEM_PATH = 'darkchem-weights/protonated'


## Functions ##

def get_latent_space_vectors(df, smiles_col_name='SMILES'):
    '''
    Adds latent space vectors to a dataframe of SMILES using DarkChem.

    Parameters
    ----------
    df : pandas dataframe
        Pandas dataframe.

    Returns
    -------
    final : pandas dataframe
        Dataframe with latent vectors.
    '''

    # save smiles to a file for darkchem use
    df[smiles_col_name].to_csv('temp_smiles.txt', index=False)

    # call darkchem from the command line
    input_ = ['darkchem predict latent temp_smiles.txt', DARKCHEM_PATH]
    input_ = ' '.join(input_)
    os.popen(input_).read()
    os.remove('temp_smiles.txt')

    # read in the file of latent space vectors
    latent_vec = np.load('temp_smiles_latent.npy')
    os.remove('temp_smiles_latent.npy')
    latent_df = pd.DataFrame(latent_vec)

    # create column names for vector dataframe
    col_names = []
    for i in range(128):
        col_names.append(('LS' + str(i + 1)))
    latent_df.columns = col_names

    # replace empty values with 'None'
    latent_df = latent_df.fillna('None')
    final = df.join(latent_df)

    # return dataframe of SMILES with latent
    # space vectors
    return final


## Main ##

if __name__ == '__main__':

    # Use library from paper as example
    fname = r'../data/SupportingData.xlsx'
    data = pd.read_excel(fname, sheet_name='SuspectLibrary',
                         usecols=['SMILES'])
    data = get_latent_space_vectors(data)
