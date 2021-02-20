
import os
import numpy as np
import pyutils.filename as fn
import pyutils.ctutils.driver as ctd
import pyutils.ctutils.flame as ctf

def counterflow_premixed_extinction(
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
    if 'a_init' in kwargs.keys():
        a_init = kwargs['a_init']
    else:
        a_init = 100.
        
    if 'a_max' in kwargs.keys():
        a_max = kwargs['a_max']
    else:
        a_max = 1.E+5

    if 'L_init' in kwargs.keys():
        L_init = kwargs['L_init']
    else:
        L_init = 0.05

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

    # iterate to get the extinction
    params['a'] = a_init

    L = L_init

    flag = 3

    while True:

        if flag == 1 :

            print('Strain rate = {:g} extinct'.format(params['a']))
            f0 /= f1
            break

        if flag == 2 :

            print('Strain rate = {:g} negative flame speed'.format(params['a']))
            #break

        if flag == 3 :

            print('Strain rate = {:g} initialization'.format(params['a']))
            solution = None

        else :

            print('Strain rate = {:g} success'.format(params['a']))
            flame_name = fn.params2name(params)
            solution = '{}.xml'.format(flame_name)

            a_old = params['a']
            L_old = L

            f0_a = np.exp( f0 )

            L = L_old / np.power(f0_a, 0.5)
            params['a'] = a_old * f0_a

            a_diff = params['a'] - a_old
            if params['a'] > a_max or a_diff < f2 * a_old :
                break

        flag = ctd.counterflow_premixed_flame(
            chemistry = chemistry,
            fuel = fuel,
            oxidizer = oxidizer,
            temperature = T,
            pressure = p,
            phi = phi,
            a = params['a'],
            solution = solution,
            width = L,
            **kwargs
        )

    os.chdir(pwd)
    
    return

def counterflow_premixed_extinction_(
    gas,
    fuel = {'CH4':1.},
    oxidizer = {'O2':1., 'N2':3.76},
    T = 300.,
    p = 101325.,
    phi = 1.,
    solution_name = 'tmp.xml',
    **kwargs):

    # read kwargs

    # parameters to approach extinction
    if 'a_init' in kwargs.keys():
        a_init = kwargs['a_init']
    else:
        a_init = 100.
        
    if 'a_max' in kwargs.keys():
        a_max = kwargs['a_max']
    else:
        a_max = 1.E+6

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

    gas.set_equivalence_ratio(phi, fuel, oxidizer)
    gas.TP = T, p
    flame = ctd.free_flame_(gas, width=L_init)

    fs_0 = ctf.PremixedFlameState(flame,fuel,oxidizer)
    T_fs = fs_0.T_peak()
    sc_0 = fs_0.consumption_speed()
    dl_0 = fs_0.thermal_thickness()

    # iterate to get the extinction
    a = a_init
    L = L_init
    solution = None

    # data containing dl, sc, I0, Ka
    d = np.array([[dl_0, sc_0],[0., 1.]])

    while True:

        gas.set_equivalence_ratio(phi, fuel, oxidizer)
        gas.TP = T, p

        flame = ctd.counterflow_premixed_flame_(
            gas,
            a = a,
            solution = solution,
            width = L,
            **kwargs
        )

        # update a and L
        f0_a = np.exp(f0)
        L /= np.power(f0_a, 0.5)
        a *= f0_a

        # solution for iteration
        flame.save(solution_name)
        solution = solution_name

        # get state
        fs_a = ctf.PremixedFlameState(flame,fuel,oxidizer,T_fs)

        sc_a = fs_a.consumption_speed()
        if sc_a < 1.0E-3:
            break

        I0_a = sc_a / sc_0
        Ka_a = fs_a.strain_rate() * dl_0 / sc_0

        #d.append([Ka_a, I0_a])
        d = np.append(d, [[Ka_a, I0_a]], axis=0)

        if a > a_max:
            break

    return d
