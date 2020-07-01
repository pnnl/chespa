# -*- coding: utf-8 -*-
'''
chespa: Streamlining chemical space evaluation of molecular sets and
assessing trends of ionizability

@author: Jamie Nunez
(C) 2020 - Pacific Northwest National Laboratory

Module for evaluating chemical formulas
'''

import re
from molmass import Formula

def is_sponch(formula):
    '''
    Checks if elements in the formula are limited to only S, P, O, N, C, and H.
    Returns True if this is the case and False if not.
    
    Parameters
    ----------
    formula: string or dict
        Molecular formula of compound. Check formula_split function desc. for limitations.
    '''
    sponch = {'S', 'P', 'O', 'N', 'C', 'H'}
    formula = formula_split(formula)
    for element in formula:
        if element not in sponch:
            return False
    return True

def match(cs_object, formula, name):
    '''
    Accepts a ChemSpider object and checks if its formula and name match the given values.
    Returns True if this is the case and False if not.
    
    Parameters
    ----------
    cs_object: object
         Object returned by ChemSpider. For example, one of the matches returned when using 
         simple_search()
    formula: string or dict
        Molecular formula of compound. Check formula_split function desc. for limitations.
    name: string
        Name of compound.
    
    Returns
    ----------
    match: boolean
        True if both the name and formula match and False if not.\n
        Note: if the ChemSpider compound doesn't have the same common name (the main name listed
        at the top of the entry page), the name check will fail and this function will return False.
    '''
    cs_formula = cs_object.molecular_formula.replace('{','').replace('}','').replace('_','')
    if type(formula) == dict: # cs_formula type must match the formula type given
        cs_formula = formula_split(cs_formula)
    formula_match = check_formula(cs_formula, formula)
    name_match = check_name(cs_object.common_name, name)
    return formula_match and name_match

def check_name(name1, name2):
    '''
    Accepts two names and returns if they are the same, accounting for the similarity between
    names like 'acetic acid' and 'acetate'. Not case sensitive. Returns True if they are
    equivalent and False otherwise.
    '''
    name1 = name1.lower()
    name2 = name2.lower()
    name1b = name1
    name2b = name2
    if name1.endswith('ate'):
        name1b = name1[:-3] + 'ic acid'
    if name2.endswith('ate'):
        name2b = name2[:-3] + 'ic acid'
    if name1 == name2 or name1b == name2 or name1 == name2b:
        return True
    else:
        return False

def check_formula(formula1, formula2):
    '''
    Accepts two formulas of the same type (string or dictionary) and checks if they are the same.
    Considers H and D to be the same element and ignores charges. Returns True if they match and
    False otherwise.
    '''
    if type(formula1) != type(formula2):
            formula1 = formula_to_string(formula1)
            formula2 = formula_to_string(formula2)
    if type(formula1) == dict:
        form1 = formula_split(formula_to_string(formula1).replace('D', 'H')) # To count D's as H's
        form2 = formula_split(formula_to_string(formula2).replace('D', 'H'))
    else:
        if '+' in formula1:
            formula1 = formula1[:formula1.index('+')]
        if '-' in formula1:
            formula1 = formula1[:formula1.index('-')]
        if '+' in formula2:
            formula2 = formula2[:formula2.index('+')]
        if '-' in formula2:
            formula2 = formula2[:formula2.index('-')]
        form1 = formula_split(formula1.replace('D', 'H')) # To count D's as H's
        form2 = formula_split(formula2.replace('D', 'H'))
    return cmp(form1, form2) == 0
        
def check_mass(mass1, mass2, max_percent_error):
    ''' 
    Accepts two mass values and returns if they are equal, within the passed percent error.
    For a percent error > 100%, use % notation, not decimal notation. Otherwise, accepts decimal
    and percent notation without error. Returns True if the masses are considered equal and False
    otherwise.
    '''
    if max_percent_error > 1:
        max_percent_error = max_percent_error / 100.0
    if mass2 * (1 + max_percent_error) > mass1 and mass2 * (1 - max_percent_error) < mass1:
        return True
    else:
        return False

def calculate_mass(formula, isotope=True):
    '''
    Accepts a formula and calculates its mass.
    
    Parameters
    ----------
    formula: string or dict
        Molecular formula of compound
    isotope: boolean
        If True, calulates the isotopic mass. If False, calculates the average mass. Default
        is True.
    '''
    formula = formula_to_string(formula)
    if isotope:
        return float(Formula(formula).isotope.mass)
    else:
        return float(Formula(formula).mass)
        
def insert_one(string):
    '''
    Helper method for formula_split. Goes through a string formula and adds ones after elements
    with an abundance of one.
    '''
    i = 0
    while i < len(string) - 1:
        if string[i].isalpha() and string[i+1].isupper():
            string = string[:i+1] + "1" + string[i+1:]
        i += 1
    if string[len(string) - 1].isalpha():
        string += "1"
    return string

def process_parenthesis(formula):
    '''
    If a formula contains parentheses, removes them. Example: C(CH3)3 as an input would return
    CCH3CH3CH3. More than one layer of parentheses causes an error.
    '''
    findparen = ''.join([i for i in formula if not i.isdigit()])
    findparen = ''.join([i for i in findparen if not i.isalpha()])
    if findparen.find('((') != -1:
        print('Error: Given formula too complex.')
        return None
    else:
        form2 = formula.split("(")
        new_formula = [form2[0]]
        for i in range(1, len(form2)):
            new = form2[i].split(")")
            for i in new:
                new_formula.append(i)
        new_formula = list(filter(None, new_formula))
        for i in range(0, len(new_formula)-1):  
            piece = str(new_formula[i+1])
            check = piece[:2]
            if check.isdigit(): # In case there is a two digit subscript after the parenthesis
                new = new_formula[i]*(int(check))
                new_formula[i] = new
                new_formula = new_formula[:i+1] + [piece[1:]] + new_formula[i+2:]
            elif check[0].isdigit():
                new = new_formula[i]*(int(check[0]))
                new_formula[i] = new
                new_formula = new_formula[:i+1] + [piece[1:]] + new_formula[i+2:]
        result = ''
        for i in new_formula:
            result += i
        return result

def process_periods(formula):
    '''
    If a formula contains periods, removes them. Example: C.3CH3 as an input would return
    CCH3CH3CH3.
    '''
    formula_list = formula.split(".")
    new_form = ""
    for piece in formula_list:
        if piece[0].isdigit():
            num = int(piece[0])
            piece = piece[1:]
            piece = piece * num
        new_form += piece
    return new_form

def _add_deut(l, formula):
    d = 0
    for piece in l:
        if 'D' in piece:
            num = piece.strip()[-1]
            if num == 'D':
                num = 1
            else:
                num = int(num)
            d += num
    formula = add_element(formula, 'D', d)
    formula = add_element(formula, 'H', -d)
    return formula
   
def inchi_formula(inchi):
    formula = formula_split(inchi.split('/')[1])
    
    # Check Isotope Layer
    isotope = [x for x in inchi.split('/') if x.startswith('i')]
    if len(isotope) > 0:
        info = isotope[0]
        if info.endswith(';'):
            info = info[:-1]
        l = info.split(',')
        if 'D' in info:
            formula = _add_deut(l, formula)
        if '+' in info:
            add = 0
            for piece in l:
                if '+' in piece:
                    temp = piece.strip()
                    index = temp.index('+')
                    add += int(temp[index:])
            formula = add_element(formula, '13C', add)
            formula = add_element(formula, 'C', -add)
    
    # Check for 2nd h layer
    isotope = [x for x in inchi.split('/') if x.startswith('h')]
    if len(isotope) > 1:
        info = isotope[1] # want the second hydrogen info
        l = info.split(',')
        if 'D' in info:
            formula = _add_deut(l, formula)
            
    # Check for missing/added protons
    isotope = [x for x in inchi.split('/') if x.startswith('p')]
    if len(isotope) > 0:
        formula = add_element(formula, 'H', int(isotope[0][1:]))
    return formula


def add_formula(dictionary, formula):
    '''
    Accepts a formula dictionary and adds a new formula to it.
    
    Parameters
    ----------
    dictionary: dict
        Formula dictionary (ex: one returned from formula_split)
    
    Returns
    ---------
    new_dictionary: dict
        Formula dictionary with the new formula added.
    '''
    new = formula_split(formula)
    for elem, abundance in new.iteritems():
        dictionary = add_element(dictionary, elem, abundance)
    return dictionary

def add_element(dictionary, element, number):
    '''
    Adds an element to passed dictionary. If the element is already there, simply adds to its
    abundance.
    
    Parameters
    ----------
    dictionary: dict
        Formula dictionary (ex: one returned from formula_split)
    element: str
        Element to be added to dictionary
    number: int
        Abundance of element to be added.
    
    Returns
    ---------
    new_dictionary: dict
        Formula dictionary with the new element added.
    '''
    if element in dictionary:
        dictionary[element] = dictionary[element] + number
    else:
        dictionary[element] = number
    if dictionary[element] == 0:
        del dictionary[element]
    return dictionary

def formula_split(formula):
    '''
    Accepts a formula in the form of a string and splits it into its elements and their
    respective abundances. Assumes the given formula has at least one element. If the passed
    formula is already in dictionary form, simply returns that same dictionary.
    
    Limitations
    ----------
    Can process formulas that contain up to one layer of parentheses and periods. More than one
    layer of parentheses, x's, and non-alphanumeric symbols will cause errors in a way that can
    cause this to crash.
    
    Parameters
    ----------
    formula: string
        Molecular formula of compound
    
    Returns
    ----------
    formula: dict
        Dictionary with key values being the elemenets present in the compound and their
        respective key values being the abundance of that element.
    '''
    if type(formula) == dict:
        return formula
    formula = str(formula)
    formula = process_parenthesis(process_periods(formula))
    if formula is not None:    
        formula = insert_one(formula)
        elements = re.split(r'[0,1,2,3,4,5,6,7,8,9]\s*', formula)
        elements = filter(None, elements)
        numbers = re.findall( r'\d+\.*\d*', formula)
        elements = list(map(str, elements))
        numbers = list(map(int, numbers))
        formula = {elements[0] : numbers[0]}
        for i in range(1, len(elements)):
            formula = add_element(formula, elements[i], numbers[i])
    return formula
    
def formula_to_string(formula):
    '''
    Accepts a formula and converts it into string form. If it is already a string, simply returns
    the same formula.
    '''
    if type(formula) == str:
        return str(formula)
    string = ''
    for key in sorted(formula):
        n = formula[key]
        if key == '13C':
            key = '(13C)'
        if n == 1:
            string += key
        else:
            string += key + str(n)
    return string