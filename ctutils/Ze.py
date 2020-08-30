# calculate the Zeldovich number by perturb nitrogen percentage in oxidizer

import os
import numpy as np
import pyutils.ctutils as pc
from pyutils.filename import params2name

def Ze(
    chemistry = 'FFCM-1.cti',
    fuel = {'CH4':1.},
    oxidizer = {'O2':1., 'N2':3.76},
    T = 300.,
    p = 1., 
    phi = 1.,
    perturb = 0.01,
    **kwargs):

    # working directory
    pwd = os.getcwd()
    
    # solution name
    flame_params = {}
    flame_params['T'] = T
    flame_params['p'] = p
    flame_params['phi'] = phi
    flame_name = params2name(flame_params)
    solution = '{}.xml'.format(flame_name)

    # perturbation
    factor = np.array([1.-perturb, 1., 1.+perturb])

    # quantaties
    flux = np.zeros(3)
    Tb = np.zeros(3)

    params = {}
    for i, f in enumerate(factor):

        # subfolder
        params['N2'] = f
        folder_name = params2name(params)

        os.makedirs(folder_name, exist_ok=True)
        os.chdir(folder_name)

        if not os.path.isfile(solution):

            # oxidizer stream with perturbation
            stream = {}
            for k, v in oxidizer.items():
                stream[k] = v
            stream['N2'] = oxidizer['N2'] * f

            pc.driver.free_flame(chemistry, fuel, stream, T, p, phi, **kwargs)

        fs = pc.flame.FreeFlameState(solution, chemistry, fuel, oxidizer)

        flux[i] = fs.mass_flux()
        Tb[i] = fs.flame.T[-1]

        os.chdir(pwd)

    grad = np.gradient(flux, Tb)
    Ze = 2.*grad[1]*(Tb[1]-T)/flux[1]

    return Ze

    