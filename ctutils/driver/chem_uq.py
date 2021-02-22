import numpy as np
import pyutils.ctutils as ctu

def free_flame(
        n = 1000,
        fuel = 'H2',
        oxidizer = {'O2':1., 'N2':3.76},
        T = 300.,
        p = 101325.,
        phi = 0.6
        ):
    
    # for unrealistic parameters
    if p < 0.:
        raise ValueError('Negative pressure')
    if T < 0.:
        raise ValueError('Negative inlet temperature')
    if phi < 0.:
        raise ValueError('Negative equivalence ratio')

    # get mech file
    mech = ctu.mechanisms.select.get_mechanism(fuel)
    # get uncertainty factors
    mech_uq = ctu.mechanisms.select.get_mech_uncertainty(mech)
    reaction_index, uncertainty_factor = ctu.sample.sample_k_factor(mech_uq,n)

    sl = np.zeros(n)
    dl = np.zeros(n)

    for i, factor in enumerate(uncertainty_factor):
        # gas object
        gas = ctu.gas.mixture(mech, fuel, oxidizer, T, p, phi)
        gas = ctu.gas.multiply(gas, reaction_index, factor)

        flame = ctu.driver.free_flame_(gas)
        fs = ctu.flame.PremixedFlameState(flame, fuel)

        sl[i] = fs.consumption_speed()
        dl[i] = fs.thermal_thickness()

    return sl, dl

def counterflow_premixed_extinction(
        n = 1000,
        fuel = {'CH4':1.},
        oxidizer = {'O2':1., 'N2':3.76},
        T = 300.,
        p = 101325.,
        phi = 1.,
        **kwargs):
    
    # for unrealistic parameters
    if p < 0.:
        raise ValueError('Negative pressure')
    if T < 0.:
        raise ValueError('Negative inlet temperature')
    if phi < 0.:
        raise ValueError('Negative equivalence ratio')

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

    # get mech file
    mech = ctu.mechanisms.select.get_mechanism(fuel)
    # get uncertainty factors
    mech_uq = ctu.mechanisms.select.get_mech_uncertainty(mech)
    reaction_index, uncertainty_factor = ctu.sample.sample_k_factor(mech_uq,n)

    d_all = []
    
    for factor in uncertainty_factor:
        # gas object
        gas = ctu.gas.mixture(mech, fuel, oxidizer, T, p, phi)
        gas = ctu.gas.multiply(gas, reaction_index, factor)

        d = ctu.driver.counterflow_premixed_extinction_(
                gas,fuel,oxidizer,T,p,phi,
                a_init=a_init,a_max=a_max,L_init=L_init)

        d_all.append(d.transpose())

    return d_all