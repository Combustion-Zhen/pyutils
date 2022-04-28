import os
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

    # IO flags
    if 'folder_overwrite' in kwargs.keys():
        flag_folder = kwargs['folder_overwrite']
    else:
        flag_folder = True

    # parameters to approach extinction
    if 'a_init' in kwargs.keys():
        a_init = kwargs['a_init']
    else:
        a_init = 100.
        
    if 'L_init' in kwargs.keys():
        L_init = kwargs['L_init']
    else:
        L_init = 0.05

    if 'L_factor' in kwargs.keys():
        L_factor = kwargs['L_factor']
    else:
        L_factor = 0.5

    # factors
    # a_{n+1} = exp(f0) * a_n
    if 'f0' in kwargs.keys():
        f0 = kwargs['f0']
    else:
        f0 = 0.1

    # f0_{n+1} = f0_n / f1
    if 'f1' in kwargs.keys():
        f1 = kwargs['f1']
    else:
        f1 = 2.0

    # a_threshold_n = a_n * f2
    if 'f2' in kwargs.keys():
        f2 = kwargs['f2']
    else:
        f2 = 1.E-3

    if 'restart' in kwargs.keys():
        rest = kwargs['restart']
    else:
        rest = False

    params = {}
    params['T'] = T
    params['p'] = p
    params['phi'] = phi

    folder_name = pu.filename.params2name(params)
    pwd = os.getcwd()

    os.makedirs(folder_name, exist_ok=flag_folder)
    os.chdir(folder_name)

    pu.ctutils.driver.free_flame(
        chemistry=chemistry, 
        fuel=fuel, oxidizer=oxidizer, 
        pressure=p, temperature=T, phi=phi,
        **kwargs)

    pressure = p * ct.one_atm
    # iterate to get the extinction
    a = a_init
    L = L_init
    a_old = 0
    L_old = L

    while True:

        gas = pu.ctutils.gas.mixture(chemistry, fuel, oxidizer, T, pressure, phi)
        flame = pu.ctutils.driver.counterflow_twin_flame(
            gas, a = a, solution = solution, width = L, **kwargs
        )

        hrr = flame.heat_release_rate.max()
        if hrr < 1.0:
            print('Strain rate {:g} extinction'.format(a))
            f0 /= f1
            f0_a = np.exp(f0)
            a = a_old * f0_a
            L = L_old / np.power(f0_a, L_factor)
            if (a - a_old) < (f2 * a_old):
                print('Iteration stop')
                break
            continue

        print('Strain rate {:g} success'.format(a))
        # solution for iteration
        params['a'] = a
        case = pu.filename.params2name(params)+'.xml'
        flame.save(case)
        solution = case
        if rest:
            solution = None

        a_old = a
        L_old = L

        # update a and L
        f0_a = np.exp(f0)
        a *= f0_a
        L /= np.power(f0_a, 0.5)

    os.chdir(pwd)
    return
