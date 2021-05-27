# opposed twin premixed flame

import numpy as np
import cantera as ct
import pyutils.ctutils.gas as cg
from pyutils.filename import params2name

def counterflow_twin_flame(
        gas,
        a = 1000.,
        solution = None,
        **kwargs
        ):

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
        ct_slope = 0.2

    if 'ct_curve' in kwargs.keys():
        ct_curve = kwargs['ct_curve']
    else:
        ct_curve = 0.2

    if 'ct_prune' in kwargs.keys():
        ct_prune = kwargs['ct_prune']
    else:
        ct_prune = 0.1

    if 'ct_max_grids' in kwargs.keys():
        ct_max_grids = kwargs['ct_max_grids']
    else:
        ct_max_grids = 5000

    # get inlet velocity based on the strain rate
    # $a_1=\dfrac{2U_1}{L}\left(1+\dfrac{U_2\sqrt{\rho_2}}{U_1\sqrt{\rho_1}}\right)$
    # $a_2=\dfrac{2U_2}{L}\left(1+\dfrac{U_1\sqrt{\rho_1}}{U_2\sqrt{\rho_2}}\right)$
    # with $\rho_1 U_1^2 = \rho_2 U_2^2$
    # $a_1=\dfrac{4U_1}{L}$ $a_2=\dfrac{4U_2}{L}$
    # set stream 1 and 2 for unburnt and equilibrium status respectively
    v = a * width / 4.

    # mass rate
    mass_flux = gas.density * v

    # Create flame object
    f = ct.CounterflowTwinPremixedFlame(gas=gas, width=width)

    f.transport_model = transport
    f.P = gas.P
    f.reactants.mdot = mass_flux

    f.set_refine_criteria(ratio=ct_ratio, 
                          slope=ct_slope, 
                          curve=ct_curve, 
                          prune=ct_prune)

    f.set_max_grid_points(f.flame, ct_max_grids)

    # load saved case if presented
    if solution is not None:

        f.restore(solution, loglevel=loglevel)

        # scaling of saved solution
        solution_width = f.grid[-1] - f.grid[0]
        width_factor = width / solution_width

        solution_a = 4.*f.u[0]/solution_width
        a_factor = a / solution_a

        normalized_grid = f.grid / solution_width

        u_factor = a_factor * width_factor

        # update solution initialization following Fiala & Sattelmayer
        f.flame.grid = normalized_grid * width
        f.set_profile('u', normalized_grid, f.u*u_factor)
        f.set_profile('V', normalized_grid, f.V*a_factor)
        f.set_profile('lambda', normalized_grid, f.L*np.square(a_factor))

        f.reactants.mdot = mass_flux

    else:

        f.set_initial_guess()

    f.solve(loglevel=loglevel, auto=True)

    return f
