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

    # parameters
    params = {}
    params['T'] = temperature
    params['p'] = pressure
    params['phi'] = phi

    case = params2name(params)

    # pressure, convert to [Pa]
    pressure *= ct.one_atm

    # gas object
    gas = cg.mixture(chemistry, fuel, oxidizer, temperature, pressure, phi)

    f = free_flame_(gas, inlet, **kwargs)

    # return for unburnt flame
    if np.max(f.T) < temperature + 100. :
        return 1

    f.save('{}.xml'.format(case))

    return 0

def free_flame_(
        gas,
        **kwargs
        ):

    # read kwargs
    if 'direct' in kwargs.keys():
        direct = kwargs['direct']
    else:
        direct = 'left'
    
    if 'transport' in kwargs.keys():
        transport = kwargs['transport']
    else:
        transport = 'Mix'

    if 'soret' in kwargs.keys():
        flag_soret = kwargs['soret']
    else:
        flag_soret = False

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
        ct_slope = 0.1

    if 'ct_curve' in kwargs.keys():
        ct_curve = kwargs['ct_curve']
    else:
        ct_curve = 0.1

    if 'ct_prune' in kwargs.keys():
        ct_prune = kwargs['ct_prune']
    else:
        ct_prune = 0.05

    if 'ct_max_grids' in kwargs.keys():
        ct_max_grids = kwargs['ct_max_grids']
    else:
        ct_max_grids = 5000

    # flame object
    f = ct.FreeFlame( gas, width=width, direct=direct )

    f.set_initial_guess()

    f.set_refine_criteria(ratio=ct_ratio, 
                          slope=ct_slope, 
                          curve=ct_curve,
                          prune=ct_prune)

    f.set_max_grid_points(f.flame, ct_max_grids)

    f.radiation_enabled = False

    try:
        f.solve( loglevel=loglevel, auto=True )
        f.transport_model = transport
        f.soret_enabled = flag_soret
        f.solve( loglevel=loglevel )
    except Exception as e:
        print('Error: not converge for case:',e)

    return f
