"""
Helper file to get substructures used in paper

Note: SPRECTRePy required and can be downloaded here:
github.com/pnnl/spectrepy

It is recommended to run this script on a supercomputer
"""

### Imports ###
import pandas as pd
import SPECTRePy


### Functions ###

def get_substructures(df):

    for ind, val in enumerate(df['Smiles']):
        Substructures = SPECTRePy.main(val)
        for j in Substructures:  # j: single substructure, will be a column

            # check if the column already exists
            if j in df.columns: # if exists, directly assigns 1 to the substructure for the cell of the corresponding Smiles
                df[j][ind] = 1  
            else:               # if not exists, creates a column filled out with 0's, then assigns 1 to the cell. 
                df[j] = [0]*len(df['Smiles'])
                df[j][ind] = 1
            
    return df 


### Main ###

if __name__ == '__main__':

    # Use library from paper as example
    fname = r'SupportingData.xlsx'
    data = pd.read_excel(fname, sheet_name='SuspectLibrary',
                         usecols=['SMILES'])
    get_substructures(data)






