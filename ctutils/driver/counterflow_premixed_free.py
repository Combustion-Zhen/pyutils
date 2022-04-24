import os
import numpy as np
import pyutils.filename as fn
import pyutils.ctutils.driver as ctd
import pyutils.ctutils.flame as ctf

def counterflow_premixed_free(
    chemistry = 'FFCM-1.cti',
    fuel = {'CH4':1.},
    oxidizer = {'O2':1., 'N2':3.76},
    T = 300.,
    p = 1.,
    phi = 1.,
    **kwargs):

    # read kwargs

    # IO flags
    if 'folder_overwrite' in kwargs.keys():
        flag_folder = kwargs['folder_overwrite']
    else:
        flag_folder = True

    # parameters to approach extinction
    if 'Ka_init' in kwargs.keys():
        Ka_init = kwargs['Ka_init']
    else:
        Ka_init = 2.
        
    if 'Ka_min' in kwargs.keys():
        Ka_min = kwargs['Ka_max']
    else:
        Ka_min = 0.02

    if 'L_init' in kwargs.keys():
        L_init = kwargs['L_init']
    else:
        L_init = 0.001

    # factors
    # a_{n+1} = exp(f0) * a_n
    if 'f0' in kwargs.keys():
        f0 = kwargs['f0']
    else:
        f0 = 0.1

    # for unrealistic parameters
    if p < 0.:
        raise ValueError('Negative pressure')
    if T < 0.:
        raise ValueError('Negative inlet temperature')
    if phi < 0.:
        raise ValueError('Negative equivalence ratio')

    params = {}
    params['T'] = T
    params['p'] = p
    params['phi'] = phi

    folder_name = fn.params2name(params)

    pwd = os.getcwd()

    os.makedirs(folder_name, exist_ok=flag_folder)
    os.chdir(folder_name)

    # calculate free flame
    ctd.free_flame(
        chemistry=chemistry, 
        fuel=fuel, oxidizer=oxidizer, 
        pressure=p, temperature=T, phi=phi,
        **kwargs)

    fs = ctf.FreeFlameState('{}.xml'.format(folder_name),
                            chemistry, fuel, oxidizer)
    sc0 = fs.consumption_speed()
    dl0 = fs.thermal_thickness()

    a = Ka_init * sc0 / dl0
    params['a'] = a
    solution = None

    L = L_init

    while True:

        ctd.counterflow_premixed_flame(
            chemistry = chemistry,
            fuel = fuel,
            oxidizer = oxidizer,
            temperature = T,
            pressure = p,
            phi = phi,
            a = a,
            solution = solution,
            width = L,
            **kwargs
        )

        flame_name = fn.params2name(params)
        solution = '{}.xml'.format(flame_name)

        f0_a = np.exp( f0 )
        L *= np.power(f0_a, 0.5)
        a /= f0_a

        if a * dl0 / sc0 < Ka_min :
            break

    os.chdir(pwd)
    
    return