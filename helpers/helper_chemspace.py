"""
Helper file to get ChemSpace parameters used in paper

Note: ChemAxon's cxcalc required, license needed.
For more info on license, visit https://chemaxon.com/discounted-licenses

Once downloaded, set CXCALC_PATH variable
"""

## Imports ##
from molmass import Formula
import pandas as pd
from subprocess import PIPE, Popen # Help with command line

from formula_module import formula_split, formula_to_string


## Globals ##
CXCALC_PATH = r'C:\Program Files\ChemAxon\MarvinSuite\bin\\'


## Functions ##

# Returns exact mass of given formula
def _get_mass(formula):

    # Potentially fix formula if too complex for Formula class 
    formula = formula_to_string(formula_split(formula))
    
    # Calculate isotopic mass and return
    return float(Formula(formula).isotope.mass)


# Returns exact masses of all formulas
# Input is df with formula_col included as a column
# Returns input df with 1 new column
def get_mass(data, formula_col='Formula', inplace=False):

    if not inplace:
        data = data.copy()

    data['Mass'] = [_get_mass(x) for x in data[formula_col]]
    return data


# Returns ratio e1/e2
def _get_ratios(formula, e1, e2):

    # Return None if formula is None
    if formula is None:
        return None
    
    # Split into dictionary of elements
    try:
        formula = formula_split(formula)
    except TypeError:
        print('Could not parse %s' % formula)
        return None

    # Get ratio
    if e2 in formula:
        if e1 in formula:
            return float(formula[e1]) / formula[e2]
        else:
            return 0
    
    # Return if e2 not in formula
    return None


# Returns the 4 element ratios used in paper for each formula
# Input is df with formula_col included as a column
# Returns input df with 4 new columns
def get_ratios(data, formula_col='Formula', inplace=False):

    if not inplace:
        data = data.copy()

    data['N/C'] = [_get_ratios(x, 'N', 'C') for x in data[formula_col]]
    data['N/H'] = [_get_ratios(x, 'N', 'H') for x in data[formula_col]]
    data['O/C'] = [_get_ratios(x, 'O', 'C') for x in data[formula_col]]
    data['O/H'] = [_get_ratios(x, 'O', 'H') for x in data[formula_col]]

    return data


# Makes call to cxcalc using function name and 2d structure
# struct_2d can be an InChI or SMILES
# Note: if you have many compounds to run, it is much faster to run a batch file
# More info about batch files: https://docs.chemaxon.com/display/docs/cxcalc+command+line+tool
def _cmdline(op, struct_2d, ref=True):
    process = Popen(
        args='cxcalc %s %s' % (op, struct_2d),
        stdout=PIPE,
        shell=True,
        cwd = CXCALC_PATH
    )
    result = str(process.communicate()[0])
    if ref:
        result = result.split('\\n')[1].split('\\t')[1].replace('\\r','')
    try:
        return float(result)
    except ValueError:
        return result


def logP(struct_2d):
    return _cmdline('logp', struct_2d)


def pka_sa(struct_2d):
    return _cmdline('pka -t acidic -a 1 -d large', struct_2d)


def bondcount(struct_2d):
    return _cmdline('bondcount', struct_2d)


def ringbondcount(struct_2d):
    return _cmdline('ringbondcount', struct_2d)


def ringbondpercent(struct_2d):
    return float(ringbondcount(struct_2d)) / bondcount(struct_2d)


def balabanindex(struct_2d):
    return _cmdline('balabanindex', struct_2d)


def hararyindex(struct_2d):
    return _cmdline('hararyindex', struct_2d)


# struct_col_name must be name of col with 2d structure (SMILES or InChI)
def get_cxcalc_properties(data, struct_col_name='SMILES', inplace=False):

    if not inplace:
        data = data.copy()

    data['Ring Bond %'] = [ringbondpercent(x) for x in data[struct_col_name]]
    data['logP'] = [logP(x) for x in data[struct_col_name]]
    data['pKa (most acidic)'] = [pka_sa(x) for x in data[struct_col_name]]
    data['Balaban Index'] = [balabanindex(x) for x in data[struct_col_name]]
    data['Harary Index'] = [hararyindex(x) for x in data[struct_col_name]]

    return

## Main ##

if __name__ == '__main__':

    # Use library from paper as example
    fname = r'../data/SupportingData.xlsx'
    data = pd.read_excel(fname, sheet_name='SuspectLibrary',
                         usecols=['Formula', 'SMILES'])
    get_mass(data, inplace=True)
    get_ratios(data, inplace=True)
    get_cxcalc_properties(data, inplace=True)  # Note: this takes significantly longer
