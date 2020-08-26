import numpy as np
import cantera as ct

class PremixedFlameState:

    def __init__(self, flame, fuel, T=None):

        self.flame = flame

        if isinstance( fuel, dict ):
            self.fuel = list(fuel.keys())
        elif isinstance( fuel, list ):
            self.fuel = fuel
        elif isinstance( fuel, str ):
            self.fuel = [fuel,]
        else:
            raise TypeError

        self.T = T

    def __idx_unburnt(self):
        return 0

    def __idx_fcr(self):
        fcr = self.fuel_consumption_rate()
        return np.argmax(fcr)

    def __idx_T(self):
        T = self.flame.T
        return np.argmax(T)

    def consumption_speed(self):

        flame = self.flame
        fuels = self.fuel

        fuel_rate = np.zeros( len(fuels) )
        fuel_mass = np.zeros( len(fuels) )

        for i, s in enumerate(fuels):

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

    def mass_flux(self):
        return self.consumption_speed()*self.flame.density[0]

    def fuel_consumption_rate(self):

        flame = self.flame
        fuels = self.fuel

        fuel_rate = np.zeros((len(fuels), flame.T.size))

        for i, s, in enumerate(fuels):

            # get species index
            index = flame.gas.species_index( s )

            fuel_rate[i] = (-flame.net_production_rates[index]
                            *flame.gas.molecular_weights[index] )

        fuel_consumption_rate = np.sum( fuel_rate, axis=0 ) 

        return fuel_consumption_rate

    def thermal_thickness(self):

        T = self.flame.T
        x = self.flame.grid

        T_grad = np.gradient( T, x )
        T_grad_max = T_grad.max()

        delta = ( T[-1] - T[0] ) / T_grad_max

        return delta

    def T_peak(self):
        return self.flame.T[self.__idx_fcr()]

    def displacement_speed(self):

        if self.T is not None:
            return np.interp(self.T, self.flame.T, self.flame.u)

        return self.flame.u[self.__idx_fcr()]

    def strain_rate(self):
        
        at = 2. * self.flame.V

        if self.T is not None:
            return np.interp(self.T, self.flame.T, at)

        return at[self.__idx_fcr()]
    
    def Ka(self):

        at = self.strain_rate()
        df = self.thermal_thickness()
        sc = self.consumption_speed()

        return at * df / sc

    def Re(self):

        df = self.thermal_thickness()
        sc = self.consumption_speed()

        rho_u = self.flame.density[0]
        mu_u = self.flame.viscosity[0]

        return sc * df * rho_u / mu_u

    def Le_eff(self, type_idx='T', type_mix='linear'):

        phi = self.equivalence_ratio()

        def mix_linear(phi):
            if phi < 0.8:
                return [1., 0.]
            elif phi > 1.2:
                return [0., 1.]
            else:
                return [3.-2.5*phi, 2.5*phi-2.]

        def mix_Bechtold(phi):
            if phi < 0.8:
                return [1., 0.]
            elif phi > 1.2:
                return [0., 1.]

        Le_spe_eff = self.Le_species_eff(type_idx)

        Le_F = self.Le_fuel(Le_spe_eff)

        Le_O = self.Le_oxidizer(Le_spe_eff)

        switch = {'linear':mix_linear, 
                  'Bechtold':mix_Bechtold}

        c = switch.get(type_mix)(phi)

        Le_eff = c[0]*Le_F+c[1]*Le_O

        return Le_eff

    def Le_fuel(self, Le_spe):

        flame = self.flame
        fuels = self.fuel

        sum_X = 0.
        sum_Le = 0.
        for i, s, in enumerate(fuels):
            # get species index
            idx = flame.gas.species_index( s )
            sum_X += flame.X[idx][0]
            sum_Le += flame.X[idx][0] * Le_spe[idx]

        Le_F = sum_Le / sum_X

        return Le_F

    def Le_oxidizer(self, Le_spe):
        
        idx = self.flame.gas.species_index('O2')

        return Le_spe[idx]

    def Le_species_eff(self, type_Le):

        Le_spe = self.Le_species()

        switch = {'T':self.__idx_T,
                  'fcr':self.__idx_fcr,
                  'unburnt':self.__idx_unburnt}

        idx = switch.get(type_Le)()

        Le_spe_eff = Le_spe[:,idx]

        return Le_spe_eff

    def Le_species(self):

        kappa = self.flame.thermal_conductivity
        cp = self.flame.cp
        rho = self.flame.density

        alpha = kappa / (rho*cp)

        D = self.flame.mix_diff_coeffs

        Le_spe = np.empty(D.shape)

        for i, D_spe in enumerate(D):
            Le_spe[i] = alpha / D_spe

        return Le_spe

    def equivalence_ratio(self):
        gas = self.flame.gas
        gas.TPY = gas.T, gas.P, self.flame.Y[:,0]
        return gas.get_equivalence_ratio()

    def export_profile(self, file_name='premix.dat', unit='cgs'):

        f = self.flame

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

class FreeFlameState(PremixedFlameState):

    def __init__(self, chemistry, fuel, solution):

        gas = ct.Solution(chemistry)
        flame = ct.FreeFlame(gas, width=0.1)

        flame.restore(solution, loglevel=0)

        PremixedFlameState.__init__(self, flame, fuel)
        
class CounterflowPremixedFlameState(PremixedFlameState):

    def __init__(self, chemistry, fuel, solution, T):

        gas = ct.Solution(chemistry)
        flame = ct.CounterflowPremixedFlame(gas, width=0.1)

        flame.restore(solution, loglevel=0)

        PremixedFlameState.__init__(self, flame, fuel, T)