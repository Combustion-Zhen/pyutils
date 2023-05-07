import numpy as np
import cantera as ct
import pyutils.ctutils as pc
import pyutils.filename as fn
from scipy.special import erfc
from scipy.signal import argrelextrema

class PremixedFlameState:

    def __init__(self, flame, fuel, oxidizer={'O2':1., 'N2':3.76}, T=None):

        self.flame = flame
        self.fuel = pc.gas.parser_stream(fuel)
        self.oxidizer = pc.gas.parser_stream(oxidizer)
        self.T = T

        if flame.T[0] < flame.T[-1]:
            self.density = flame.density[0]
        else:
            self.density = -flame.density[-1]

    def __idx_unburnt(self):
        return 0

    def __idx_fcr(self):
        fcr = self.fuel_consumption_rate()
        return np.argmax(fcr)

    def __idx_hrr(self):
        hrr = self.flame.heat_release_rate
        return np.argmax(hrr)

    def __idx_T(self):
        T = self.flame.T
        return np.argmax(T)

    def idx_ref(self):
        flag = 0
        u = self.flame.velocity
        # find the local minima
        idx = argrelextrema(u, np.less)
        if not idx[0].any():
            flag = 1
            x = self.flame.grid
            du = np.gradient(u, x)
            idx = argrelextrema(du, np.greater)
            if not idx[0].any():
                return 0, 2
        return idx[0][0], flag

    def fuel_list(self):
        return list(self.fuel.keys())

    def expansion(self):
        return self.flame.density[-1]/self.flame.density[0]

    def sigma_hrr(self):
        return self.flame.density[0]/self.flame.density[self.__idx_hrr()]

    def sigma_fcr(self):
        return self.flame.density[0]/self.flame.density[self.__idx_fcr()]

    def consumption_speed(self):

        flame = self.flame
        fuels = self.fuel_list()

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

        sc = fuel_rate_sum / ( self.density * fuel_mass_sum )

        return sc

    def thermal_consumption_speed(self):

        flame = self.flame

        T_u = flame.T[0]
        T_b = flame.T[-1]

        dT = flame.heat_release_rate/flame.cp
        sum = np.trapz(dT, flame.grid)

        sc = sum / (T_b-T_u) / self.density

        return sc

    def oxidizer_consumption_speed(self):

        flame = self.flame

        index = flame.gas.species_index('O2')

        rate = - ( np.trapz(flame.net_production_rates[index], flame.grid)
                  *flame.gas.molecular_weights[index] )

        mass = flame.Y[index, 0] - flame.Y[index,-1]

        sc = rate / (self.density * mass)

        return sc

    def mass_flux(self):
        return self.consumption_speed()*self.flame.density[0]

    def fuel_consumption_rate(self):

        flame = self.flame
        fuels = self.fuel_list()

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

        if T[-1] > T[0]:
            return ( T[-1] - T[0] ) / T_grad.max()
        else:
            return ( T[-1] - T[0] ) / T_grad.min()

    def diffusive_thickness(self):

        kappa = self.flame.thermal_conductivity
        cp = self.flame.cp
        rho = self.flame.density

        alpha = kappa / (rho*cp)

        sc = self.consumption_speed()

        if self.flame.T[0] < self.flame.T[-1] :
            return alpha[0] / sc
        else:
            return alpha[-1] / sc

    def T_peak(self):
        return self.flame.T[self.__idx_fcr()]

    def flame_location(self, T=None, direction='inward'):

        if T is None and self.T is None:
            return self.flame.grid[self.__idx_hrr()]
        elif T is not None:
            T_flame = T
        else:
            T_flame = self.T

        if direction == 'inward':
            flameT = self.flame.T
            flameG = self.flame.grid
        else:
            flameT = self.flame.T[::-1]
            flameG = self.flame.grid[::-1]

        return np.interp(T_flame, flameT, flameG)


    def displacement_speed(self, T=None, direction='inward'):

        if T is None and self.T is None:
            return self.flame.velocity[self.__idx_hrr()]
        elif T is not None:
            T_flame = T
        else:
            T_flame = self.T

        if direction == 'inward':
            flameT = self.flame.T
            flameV = self.flame.velocity
        else:
            flameT = self.flame.T[::-1]
            flameV = self.flame.velocity[::-1]

        return np.interp(T_flame, flameT, flameV)

    def displacement_speed_ref(self):
        idx, flag = self.idx_ref()
        return self.flame.velocity[idx], flag

    def density_weighted_displacement_speed(self, T=None):

        if T is not None:
            x = np.interp(T, self.flame.T, self.flame.grid)
        elif self.T is not None:
            x = np.interp(self.T, self.flame.T, self.flame.grid)
        else:
            x = self.flame.grid[self.__idx_hrr()]

        rho = np.interp(x, self.flame.grid, self.flame.density)
        sd = np.interp(x, self.flame.grid, self.flame.velocity)
            
        return rho*sd/self.density

    def strain_rate(self, T=None):
        
        if float(ct.__version__[:3]) <=2.4:
            at = 2. * self.flame.V
        else:
            at = 2. * self.flame.spread_rate

        if T is not None:
            return np.interp(T, self.flame.T, at)
        elif self.T is not None:
            return np.interp(self.T, self.flame.T, at)
        else:
            return at[self.__idx_hrr()]

    def strain_rate_ref(self):
        idx, flag = self.idx_ref()
        if flag == 2:
            return 0
        return 2.*self.flame.spread_rate[idx]
    
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

    def Le_fuel(self, Le_spe):

        flame = self.flame
        fuels = self.fuel_list()

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
        return gas.equivalence_ratio()

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

    def export_profile_YT(self, file_name='cema.inp'):
        f = self.flame
        data = np.zeros((f.grid.size, len(f.gas.species_names)+1))
        data[:,:-1] = f.Y.transpose()
        data[:,-1] = f.T
        np.savetxt(file_name, data, fmt='%20.10E', delimiter='')
        return 0

class FreeFlameState(PremixedFlameState):

    def __init__(self, solution, chemistry, fuel, oxidizer={'O2':1., 'N2':3.76}):

        self.chemistry = chemistry

        gas = ct.Solution(chemistry, loglevel=0)
        flame = ct.FreeFlame(gas, width=0.1)

        flame.restore(solution, loglevel=0)

        PremixedFlameState.__init__(self, flame, fuel, oxidizer)

    """
    def Ze(self, perturb=0.01, **kwargs):

        chemistry = self.chemistry
        fuel = self.fuel
        oxidizer = self.oxidizer
        T = self.flame.T[0]
        p = self.flame.P / ct.one_atm
        phi = self.equivalence_ratio()

        return pc.Ze(chemistry, fuel, oxidizer, T, p, phi, perturb, **kwargs)
    """
    def Ze(self):

        tempu = self.flame.T[0]
        tempb = self.flame.T[-1]
        ddT = np.gradient(np.gradient(self.flame.T, self.flame.grid), self.flame.grid)
        temp0 = self.flame.T[ddT.argmax()]

        return 4*(tempb-tempu)/(tempb-temp0)

    def nTb(self):
        
        tempu = self.flame.T[0]
        tempb = self.flame.T[-1]

        return tempb / ( tempb - tempu)

    """
    calculate the difference of Ma between the consumption speed and displacement speed
    Giannakopoulos et al., CNF, 2019
    """
    def Ma_diff(self):

        T = self.flame.T
        kappa = self.flame.thermal_conductivity

        ns = T[-1] / T[0]
        nT = T / T[0]
        nL = kappa / kappa[0]

        fcr = self.fuel_consumption_rate()
        idx = np.argmax(fcr)

        func_0 = nL/nT
        func_1 = nL/(nT-1)

        int_0 = np.trapz(func_0, nT) * ns / (ns-1)
        int_1 = np.trapz(func_0[:idx+1], nT[:idx+1])
        int_2 = np.trapz(func_1[idx:], nT[idx:])

        L_diff_curv = int_0 - int_1
        L_diff_strain = L_diff_curv - int_2

        Ma_diff_curv = L_diff_curv * self.diffusive_thickness() / self.thermal_thickness()
        Ma_diff_strain = L_diff_strain * self.diffusive_thickness() / self.thermal_thickness()

        return Ma_diff_strain, Ma_diff_curv

    def Ma(self):
        
        return 1/self.Le_eff() + self.Ze()/2*(1-1/self.Le_eff())
        
    def Le_eff(self, type_idx='unburnt', type_mix='erf'):

        def mix_linear(Le_F, Le_O, phi):
            if phi < 0.8:
                return Le_F
            elif phi > 1.2:
                return Le_O
            else:
                return (3.-2.5*phi)*Le_F+(2.5*phi-2.)*Le_O

        def mix_erf(Le_F, Le_O, phi):
            Ze = self.Ze()
            phi_n = phi/(1.+phi)
            x = Ze*(phi_n-0.5)*2.
            f = erfc(x)
            return Le_O + (Le_F-Le_O)*f/2.

        def mix_Bechtold(Le_F, Le_O, phi):

            Ze = self.Ze()

            if phi < 1.:
                phi_ = 1./phi
                Le_E = Le_O
                Le_D = Le_F
            else:
                phi_ = phi
                Le_E = Le_F
                Le_D = Le_O

            A = 1. + Ze * ( phi_ - 1. )

            return (Le_E+Le_D*A)/(1.+A)

        def mix_Bechtold_cut(Le_F, Le_O, phi):
            
            if phi < 0.8:
                return Le_F
            elif phi > 1.2:
                return Le_O
            else:
                return mix_Bechtold(Le_F, Le_O, phi)

        def mix_Dortz(Le_F, Le_O, phi):
            
            if phi <= 0.6:
                return Le_F
            elif phi >= 1.2:
                return Le_O
            else:
                Le_BM = mix_Bechtold(Le_F, Le_O, phi)
                if phi <= 1.:
                    return 2.5*(1.-phi)*Le_F+(2.5*phi-1.5)*Le_BM
                else:
                    return 2.5*(phi-1.)*Le_O+(3.5-2.5*phi)*Le_BM

        phi = self.equivalence_ratio()

        Le_spe_eff = self.Le_species_eff(type_idx)

        Le_F = self.Le_fuel(Le_spe_eff)

        Le_O = self.Le_oxidizer(Le_spe_eff)

        switch = {'linear':mix_linear, 
                'erf':mix_erf,
                'Bechtold':mix_Bechtold,
                'Bechtold_cut':mix_Bechtold_cut,
                'Dortz':mix_Dortz}

        Le_eff = switch.get(type_mix)(Le_F, Le_O, phi)

        return Le_eff

class ForcedPolarFlameState(PremixedFlameState):
    
    def __init__(self, solution, chemistry, fuel, oxidizer={'O2':0.21, 'N2':0.79}, direct='inward', T=None):
        
        self.chemistry = chemistry

        gas = ct.Solution(chemistry, loglevel=0)
        flame = ct.ForcedPolarFlame(gas, width=0.2, direct=direct)

        flame.restore(solution, loglevel=0)

        PremixedFlameState.__init__(self, flame, fuel, oxidizer, T)


class FreePolarFlameState(PremixedFlameState):
    
    def __init__(self, solution, chemistry, fuel, oxidizer={'O2':0.21, 'N2':0.79}, T=None):
        
        self.chemistry = chemistry

        gas = ct.Solution(chemistry, loglevel=0)
        flame = ct.FreePolarFlame(gas, width=0.2)

        flame.restore(solution, loglevel=0)

        PremixedFlameState.__init__(self, flame, fuel, oxidizer, T)

class CounterflowPremixedFlameState(PremixedFlameState):

    def __init__(self, solution, chemistry, fuel, oxidizer={'O2':1., 'N2':3.76}, T=None):

        self.chemistry = chemistry

        gas = ct.Solution(chemistry, loglevel=0)
        flame = ct.CounterflowPremixedFlame(gas, width=0.1)

        flame.restore(solution, loglevel=0)

        PremixedFlameState.__init__(self, flame, fuel, oxidizer, T)

class CounterflowTwinFlameState(PremixedFlameState):

    def __init__(self, solution, chemistry, fuel, oxidizer={'O2':1., 'N2':3.76}, T=None):

        self.chemistry = chemistry

        gas = ct.Solution(chemistry, loglevel=0)
        flame = ct.CounterflowTwinPremixedFlame(gas, width=0.1)

        flame.restore(solution, loglevel=0)

        PremixedFlameState.__init__(self, flame, fuel, oxidizer, T)
