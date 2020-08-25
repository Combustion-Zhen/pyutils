import cantera as ct
import numpy as np
import pyutils.ctutils.gas as cg
from pyutils.filename import params2name

def free_flame(
        chemistry = 'gri30.cti',
        fuel = {'CH4':1.},
        oxidizer = {'O2':1., 'N2':3.76},
        temperature = 300.,
        pressure = 1.,
        phi = 1.,
        **kwargs
        ):

    # for unrealistic parameters
    if pressure < 0.:
        raise ValueError('Negative pressure')
    if temperature < 0.:
        raise ValueError('Negative inlet temperature')
    if phi < 0.:
        raise ValueError('Negative equivalence ratio')

    # read kwargs
    if 'transport' in kwargs.keys():
        transport = kwargs['transport']
    else:
        transport = 'Mix'

    if 'width' in kwargs.keys():
        width = kwargs['width']
    else:
        width = 0.05

    if 'loglevel' in kwargs.keys():
        loglevel = kwargs['loglevel']
    else:
        # supress log output
        loglevel = 0

    # kwargs for flame solver
    if 'ct_ratio' in kwargs.keys():
        ct_ratio = kwargs['ct_ratio']
    else:
        ct_ratio = 2.

    if 'ct_slope' in kwargs.keys():
        ct_slope = kwargs['ct_slope']
    else:
        ct_slope = 0.02

    if 'ct_curve' in kwargs.keys():
        ct_curve = kwargs['ct_curve']
    else:
        ct_curve = 0.02

    if 'ct_prune' in kwargs.keys():
        ct_prune = kwargs['ct_prune']
    else:
        ct_prune = 0.01

    # parameters
    params = {}
    params['T'] = temperature
    params['p'] = pressure
    params['phi'] = phi

    case = params2name(params)

    # pressure, convert to [Pa]
    pressure *= ct.one_atm

    # gas object
    #gas = ct.Solution(chemistry)
    # construct mixture
    #mixture = cg.mixture_two_streams( gas, fuel, oxidizer, phi )
    # assign inlet gas properties
    #gas.TPX = temperature, pressure, mixture
    gas = cg.mixture(chemistry, fuel, oxidizer, temperature, pressure, phi)

    # flame object
    f = ct.FreeFlame( gas, width=width )

    f.set_initial_guess(locs=np.linspace(0., 1., num=10))

    f.set_refine_criteria(ratio=ct_ratio, 
                          slope=ct_slope, 
                          curve=ct_curve,
                          prune=ct_prune)

    f.soret_enabled = False
    f.radiation_enabled = False

    f.transport_model = transport

    try:
        f.solve( loglevel=loglevel, auto=True )
    except Exception as e:
        print('Error: not converge for case:',e)
        return -1

    # return for unburnt flame
    if np.max(f.T) < temperature + 100. :
        return 1

    f.save('{}.xml'.format(case))

    return 0

