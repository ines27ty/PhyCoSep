from math import *
from pylab import *
from ct import *
#matplotlib inline 

I = 1.0e-2                        # mol/L
a = 100.0e-9                      # m
A = 1.0e-19                       # J
eta   = 0.001                     # Pa.s
n     = 4.45e15                   # /m^3
psi_s = 0.050 * elc/kT            # dimensionless potential  psi = Psi * e / kT
z     = 1.0                       # monovalent ions 

for psi_s in (-50.0e-3* elc/kT , -25.0e-3* elc/kT): 
    for I in (1.0e-4, 1.0e-3, 1.0e-2):

        print("\n--------------------------------------------")

        print("I=", I, "psi_s=", psi_s, "\n---------------------------")
        n0 = I*1000*Na
        print("n0 = ", n0, "molecules/m^3")
        lD = sqrt(eps*kT/(2*n0*elc**2))
        print("Screening length = ", lD)

        kappa = 1/lD
        print("kappa = ", kappa)
        ka = kappa*a
        print("ka = ", ka)

        # Derjaguin approximation for the electrostatic and vdW forces, valid if ka>>1 and h<<a
        F0_EL = eps*(kT/elc)**2*32*pi*ka*tanh(psi_s/4)**2
        F0_vdW = -A*a/12

        # Load the non-linear equation solver "fsolve"
        from scipy.optimize import fsolve
        # solve "sum of forces = 0"
        res = fsolve( lambda h: F0_EL*exp(-kappa*h) + F0_vdW/h**2, 2.0e-9 )
        # fsolve returns an array with 1 value. The following line is just to extract this value
        h = res[0]
        print("h =", h )

        # Define the potentials at the location of the max h 
        V_EL_max = eps*(kT/elc)**2*32*pi*a*tanh(psi_s/4)**2*exp(-kappa*h)
        V_vdW_max = -A*a/(12*h)
        Vmax = V_EL_max + V_vdW_max
        print("V max = ", Vmax/kT)

        # Compute the stability ratio 
        W = 1 / (2*kappa*a) * exp(Vmax/kT)
        print("W     =", W)

        # perikinetic aggregation constant for perfectly effective interactions 
        k = 8*kT/(3*eta)

        # real aggregation rate constant 
        if W>1:
            kstar = k/W 
        else:
            kstar = k

        # aggregation time 
        print("tau    =", 2/(kstar*n), 's' )
        
        # CCC from formula 
        gamma0 = tanh(psi_s/4.) 
        print(gamma0 )
        CCC = 3.18e-36*gamma0**4/A**2/z**6 
        print("CCC =", CCC/1000, "mol/L") 