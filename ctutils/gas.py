import numpy as np
import cantera as ct
from pyutils.filename import name2params

def mixture(chemistry, fuel, oxidizer, T, P, phi):

    gas = ct.Solution(chemistry)
    gas.set_equivalence_ratio(phi, fuel, oxidizer)
    gas.TP = T, P

    return gas

def multiply(gas, reaction_index, uncertainty_factor):
    for k, v in enumerate(reaction_index):
        gas.set_multiplier(uncertainty_factor[k], v)
    return gas

def parser_stream(stream):

    if isinstance( stream, dict ):
        return stream
    elif isinstance( stream, str ):
        return name2params(stream, s1=',', s2=':', default=1.)

def mixture_two_streams(gas, fuel, oxidizer, phi):
    # get mixture by two streams and equivalence ratio

    # get stoichiometric coefficient

    # fuel stream
    nu_fuel = stoichiometric_nu( gas, fuel )
    if nu_fuel <= 0.:
        raise ValueError('Fuel stream')

    # oxidizer stream
    nu_oxidizer = stoichiometric_nu( gas, oxidizer )
    if nu_oxidizer >= 0.:
        raise ValueError('Oxidizer stream')

    nu = -nu_fuel/nu_oxidizer

    # construct mixture

    mixture = {}
    # add fuel
    for k, v in fuel.items():
        mixture[k] = v * phi
    # add oxidizer
    for k, v in oxidizer.items():
        if k in mixture.keys():
            mixture[k] += v * nu
        else:
            mixture[k] = v * nu

    return mixture

def stoichiometric_nu(gas, stream):
    # get the stoichiometric coefficient for a stream
    # input:
    #   gas:    cantera solution object
    #   stream: dictionary containing species name and mole number
    #           {'CH4':1, 'O2':1, 'N2':3.76}

    atom_list = ['C', 'H', 'O']
    atom_rate = [1, 0.25, -0.5]

    nu_stream = 0.

    for k, v in stream.items():

        index = gas.species_index(k)

        nu = 0.
        for i, atom in enumerate(atom_list):
            if atom in gas.element_names:
                nu += gas.n_atoms(index, atom) * atom_rate[i]

        nu_stream += nu * v

    return nu_stream

def specific_mole_number_element(gas, element):
    """
    calculate the specific mole number of an element in a mixture
    """

    return gas.elemental_mass_fraction(element)/gas.atomic_weight(element)

def specific_mole_number_elements(gas, definition='Bilger'):
    """
    calculate the specific mole number of elements with 
    the coefficient by Bilger 2 C 0.5 H -1 O
    """

    if definition == 'Bilger':
        elements = 'CHO'
    else:
        elements = definition

    z = 0.0
    for element in elements:
        if element == 'C':
            c = 2.0
        elif element == 'H':
            c = 0.5
        elif element == 'O':
            c = -1.0
        else:
            raise ValueError('Invalid element name')
        z += c*specific_mole_number_element(gas, element)

    return z

def mixture_fraction(gas, fuel, oxidizer, definition='Bilger'):

    zG = specific_mole_number_elements(gas, definition)
    zF = specific_mole_number_elements(fuel, definition)
    zO = specific_mole_number_elements(oxidizer, definition)

    return (zG-zO)/(zF-zO)

def progress_variable_normalized(gas, unburnt, burnt, definition='4spe'):

    cG = progress_variable(gas, definition)
    cU = progress_variable(unburnt, definition)
    cB = progress_variable(burnt, definition)

    return (cG-cB)/(cU-cB)

def progress_variable(gas, definition=['CO2', 'CO', 'H2O', 'H2']):

    if definition == 'T':
        return gas.T
    else:
        c = 0.0
        for spe in definition:
            c += gas.Y[gas.species_index(spe)]
        return c
