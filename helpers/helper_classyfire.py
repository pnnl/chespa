# -*- coding: utf-8 -*-
"""
Helper file to get ClassyFire superclass, as used in paper

Calls ClassyFire website and pulls Superclass info
"""

## Imports ##
from bs4 import BeautifulSoup
import pandas as pd
import requests
from time import sleep


## Globals ##
CB_URL = 'http://classyfire.wishartlab.com/entities/%s'


## Functions ##

# Calls ClassyFire site using InChIKey.
# Does not submit new compounds (must have been submitted in the past)
# Level can be any of those defined by ClassyFire
def classyfire_query(inchi_key, level='Superclass'):
    sleep(0.25)  # To reduce strain on site

    try:
        response = requests.get(CB_URL % inchi_key)
    except requests.ConnectionError:
        sleep(600)  # Give the site a 10 min break
        response = requests.get(CB_URL % inchi_key)
        
    soup = BeautifulSoup(response.content)
    result = [a.text.strip() for a in soup.select('li')]

    if len(result) > 0:
        info = [x for x in result if x.startswith(level + '\n')]
        if len(info) > 0:
            info = [x.strip() for x in info[0].split('\n') if len(x) > 0][1]
            return info
    return None


# Queries ClassyFire site for Superclass of given InChIKeys
# Does not submit new compounds (must have been submitted in the past)
# Input is df with formula_col included as a column
# Returns input df with 1 new column
def get_superclass(data, struct_col_name='InChIKey', inplace=False):

    if not inplace:
        data = data.copy()

    data['SuperClass'] = [classyfire_query(x) for x in data[struct_col_name]]
    return data


## Main ##

if __name__ == '__main__':
    
    # Use library from paper as example
    fname = r'../data/SupportingData_v5.xlsx'
    data = pd.read_excel(fname, sheet_name='SuspectLibrary', usecols=['InChIKey'])
    get_superclass(data, inplace=True)
