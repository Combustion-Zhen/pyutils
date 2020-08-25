#%%
import numpy as np
import cantera as ct
import pyutils.ctutils.gas as cg
from pyutils.filename import params2name

# %%

# one stream unburnt, one stream equilibrium

def counterflow_premixed_flame(
        chemistry = 'gri30.xml',
        fuel = {'CH4':1.},
        oxidizer = {'O2':1, 'N2':3.76},
        temperature = 300.,
        pressure = 1.,
        phi = 1.,
        a = 1000.,
        solution = None,
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

    if 'ct_max_grids' in kwargs.keys():
        ct_max_grids = kwargs['ct_max_grids']
    else:
        ct_max_grids = 5000

    # case name
    params = {}
    params['T'] = temperature
    params['p'] = pressure
    params['phi'] = phi
    params['a'] = a

    case = params2name(params)

    # pressure
    pressure *= ct.one_atm

    # gas object
    #gas = ct.Solution(chemistry)
    # construct mixutre
    #mixture = cg.mixture_two_streams(gas, fuel, oxidizer, phi)
    # unburnt stream
    #gas.TPX = temperature, pressure, mixture
    gas = cg.mixture(chemistry, fuel, oxidizer, temperature, pressure, phi)

    rho_u = gas.density

    # burnt stream
    gas.equilibrate('HP')
    rho_b = gas.density

    gas = cg.mixture(chemistry, fuel, oxidizer, temperature, pressure, phi)

    # get inlet velocity based on the strain rate
    # $a_1=\dfrac{2U_1}{L}\left(1+\dfrac{U_2\sqrt{\rho_2}}{U_1\sqrt{\rho_1}}\right)$
    # $a_2=\dfrac{2U_2}{L}\left(1+\dfrac{U_1\sqrt{\rho_1}}{U_2\sqrt{\rho_2}}\right)$
    # with $\rho_1 U_1^2 = \rho_2 U_2^2$
    # $a_1=\dfrac{4U_1}{L}$ $a_2=\dfrac{4U_2}{L}$
    # set stream 1 and 2 for unburnt and equilibrium status respectively
    v_u = a * width / 4.0
    v_b = np.sqrt( rho_u*np.square(v_u) / rho_b )

    # mass rate
    m_u = rho_u * v_u
    m_b = rho_b * v_b

    # Create flame object
    f = ct.CounterflowPremixedFlame(gas=gas, width=width)

    f.transport_model = transport
    f.P = pressure
    f.reactants.mdot = m_u
    f.products.mdot = m_b

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

        f.reactants.mdot = m_u
        f.products.mdot = m_b

    else:

        f.set_initial_guess()

    f.solve(loglevel=loglevel, auto=True)

    HRR = f.heat_release_rate

    idx = HRR.argmax()

    if HRR[idx] > 1000 :

        f.save('{}.xml'.format(case))

        if f.u[idx] > 0 :

            return 0

        else :

            return 2

    else:

        return 1