from numpy import pi 
# Boltzmann const x Temperature
kT  = 1.3806503e-23*293
# electron charge
elc = 1.60217646e-19
# Avogadro number
Na  = 6.0221415e23
# eps0
eps0 = 8.854187817620e-12
# Eps_r x Eps_0 for water
eps = 78*eps0
# Bjerum length
lambdaB = elc**2/(4.0*pi*kT*eps)
lB = lambdaB 
