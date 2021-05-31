import numpy as np
import cantera as ct
import pyutils as pu

def counterflow_twin_extinction(
    chemistry = 'FFCM-1.cti',
    fuel = {'CH4':1.},
    oxidizer = {'O2':1., 'N2':3.76},
    T = 300.,
    p = 1.,
    phi = 1.,
    solution = None,
    **kwargs):

    # read kwargs

    # parameters to approach extinction
    if 'a_init' in kwargs.keys():
        a_init = kwargs['a_init']
    else:
        a_init = 100.
        
    if 'L_init' in kwargs.keys():
        L_init = kwargs['L_init']
    else:
        L_init = 0.05

    # factors
    # a_{n+1} = exp(f0) * a_n
    if 'f0' in kwargs.keys():
        f0 = kwargs['f0']
    else:
        f0 = 0.2

    params = {}
    params['T'] = T
    params['p'] = p
    params['phi'] = phi

    pressure = p * ct.one_atm

    gas = pu.ctutils.gas.mixture(chemistry, fuel, oxidizer, T, pressure, phi)

    flame = pu.ctutils.driver.free_flame_(gas, width=L_init)

    case = pu.filename.params2name(params)+'.xml'

    flame.save(case)

    # iterate to get the extinction
    a = a_init
    L = L_init

    while True:

        gas = pu.ctutils.gas.mixture(chemistry, fuel, oxidizer, T, pressure, phi)

        flame = pu.ctutils.driver.counterflow_twin_flame(
            gas,
            a = a,
            solution = solution,
            width = L,
            **kwargs
        )

        hrr = flame.heat_release_rate.max()
        if hrr < 1.0:
            break

        params['a'] = a
        case = pu.filename.params2name(params)+'.xml'

        # solution for iteration
        flame.save(case)
        solution = case

        # update a and L
        f0_a = np.exp(f0)
        L /= np.power(f0_a, 0.5)
        a *= f0_a

    return
