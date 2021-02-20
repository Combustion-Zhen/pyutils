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
