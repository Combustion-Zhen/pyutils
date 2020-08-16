import numpy as np
import cantera as ct

class PremixedFlameState:

    def __init__(self, flame, fuel, T=None):

        sc, df, T, sd, K = premixed_flame_state(flame, fuel, T=T)

        self.sc = sc
        self.df = df
        self.T_peak = T
        self.sd = sd
        self.K = K

def premixed_flame_state(flame, fuel, T=None):
    
    sc = consumption_speed(flame, fuel)
    df = thermal_thickness(flame)

    fcr = fuel_consumption_rate(flame, fuel)
    idx = np.argmax(fcr)

    T_peak = flame.T[idx]
    sd = flame.u[idx]

    strain_rate = 2. * flame.V
    K = strain_rate[idx]

    if T is not None:
        sd = np.interp(T, flame.T, flame.u)
        K = np.interp(T, flame.T, strain_rate)

    return sc, df, T_peak, sd, K

def consumption_speed(flame, fuel):

    # check the fuel info
    if isinstance( fuel, str ):
        # single component
        fuel_list = [fuel,]
    elif isinstance( fuel, list ):
        fuel_list = fuel

    fuel_rate = np.zeros( len(fuel_list) )
    fuel_mass = np.zeros( len(fuel_list) )

    for i, s in enumerate(fuel_list):

        # get species index
        index = flame.gas.species_index( s )

        # calculate fuel consumption
        fuel_rate[i] = - ( np.trapz(flame.net_production_rates[index],
                                    flame.grid)
                          *flame.gas.molecular_weights[index] )

        # fuel mass fraction difference
        fuel_mass[i] = flame.Y[index, 0] - flame.Y[index,-1]

    fuel_rate_sum = np.sum( fuel_rate )
    fuel_mass_sum = np.sum( fuel_mass )

    sc = fuel_rate_sum / ( flame.density[0] * fuel_mass_sum )

    return sc

def fuel_consumption_rate(flame, fuel):

    # check the fuel info
    if isinstance( fuel, str ):
        # single component
        fuel_list = [fuel,]
    elif isinstance( fuel, list ):
        fuel_list = fuel

    fuel_rate = np.zeros((len(fuel_list), flame.T.size))

    for i, s, in enumerate(fuel_list):

        # get species index
        index = flame.gas.species_index( s )

        fuel_rate[i] = (-flame.net_production_rates[index]
                        *flame.gas.molecular_weights[index] )

    fuel_consumption_rate = np.sum( fuel_rate, axis=0 ) 

    return fuel_consumption_rate

def thermal_thickness(flame):

    T_grad = np.gradient( flame.T, flame.grid )
    T_grad_max = T_grad.max()

    delta = ( flame.T[-1] - flame.T[0] ) / T_grad_max

    return delta

def export_profile(f, file_name='premix.dat', unit='cgs'):

    # export cantera flame result in the form of premix output
    # variable names    (A20)
    # X, U, RHO, Y, T   (E20.10)

    # unit system
    # SI
    convertor_length = 1.0
    convertor_density = 1.0
    # cgs
    if unit == 'cgs':
        convertor_length = 1.0E+02
        convertor_density = 1.0E-03

    # variale names
    species_names =  f.gas.species_names
    variable_names = ['X', 'U', 'RHO'] + species_names +['T',]
    str_names = ''.join(['{:>20}'.format(n) for n in variable_names])

    # data for output
    data = np.zeros((f.grid.size, len(variable_names)))

    data[:,0] = f.grid * convertor_length
    data[:,1] = f.u * convertor_length
    data[:,2] = f.density * convertor_density
    data[:,3:-1] = f.Y.transpose()
    data[:,-1] = f.T

    np.savetxt(file_name, data, fmt='%20.10E', delimiter='', 
               header=str_names, comments='')

    return 0
