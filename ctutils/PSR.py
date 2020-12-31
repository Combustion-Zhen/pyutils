# %%
import cantera as ct
# %%

gas = ct.Solution('h2_Burke_n2.cti')

r = ct.IdealGasConstPressureReactor(gas)
# %%

dt = 1e-5

t = 0.

# %%
r.T
t += dt
sim = ct.ReactorNet([r])
sim.advance(t)
r.T
# %%

gas.TPX = 300, 101325, {'H':1.}
r.syncState()
r.T
# %%
sim.reinitialize()
# %%
t += dt
sim.advance(t)
r.T
# %%