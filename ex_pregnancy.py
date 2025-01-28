from math import *
from pylab import *
#from ct import *            # File of constants
import numpy as np  
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constants
eps = 6.90e-10      # F/m
elc = 1.60e-19      # C
Na = 6.02e23        # mol^-1
k = 1.38e-23        # J/K
T = 273.15+25       # K
a = 100e-9          # m
n = 4.45e15         # m^-3
A = 1.0e-19         # J
eta   = 0.001                     # Pa.s

# Different cases
I_list = [1e-4 , 1e-3, 1e-2]                # M (mol/L)
Vmax_absence = [153.9,100.12,17.6]          # Vmax  pour zeta = -50 mV
Vmax_presence = [26, 2.1, -2.1421679454330065e-47]                 # Vmax  pour zeta = -20 mV
W_50_list = []
W_25_list = []
psi_s_absence = -50e-3*elc/(k*T)
psi_s_presence = -25e-3*elc/(k*T)
z = 1                                       # Charge of the ions
k_peri = 8*k*T/(3*eta)                           # perikinetic aggregation constant for perfectly effective interactions 

for i in range(len(I_list)) :
    ## Absence of pregnancy
    print("\n--------------------------------------------")
    print("I (M = mol/L)=", I_list[i], "and Vmax_50 = ", Vmax_absence[i])
    # Calculate the Debye length and kappa
    lambda_d = np.sqrt(eps*k*T/(2*1000*Na*I_list[i]*elc**2)) 
    print("lambda_d(m) = ", lambda_d)
    kappa = 1/lambda_d
    print("kappa = ", kappa)

    # Calculate the stability ratio
    W_50 = (1/(kappa*2*a))*exp(Vmax_absence[i])    
    print("W_50 (absence) = ", W_50)
    W_50_list.append(W_50)

    # Ionic strength
    F_EL_0_a = eps * (k*T/elc)**2 * 32 * pi * kappa*a * tanh(psi_s_absence/4)**2
    F_vdw_0_a = - A * a /12
    res = fsolve( lambda h_a: F_EL_0_a*exp(-kappa*h_a) + F_vdw_0_a/h_a**2, 2.0e-9 )
    # fsolve returns an array with 1 value. The following line is just to extract this value
    h_a = res[0]
    print("h =", h_a )
    
    # CCC absence
    gamma_0_absence = tanh(psi_s_absence/4.)
    print("gamma_0_absence = ", gamma_0_absence)
    CCC_absence = 3.8e-36*(gamma_0_absence**4)/(A**2*z**6)
    print("CCC absence (mol/m3) = ", CCC_absence)
    
     # Define the potentials at the location of the max h 
    V_EL_max_a = eps*(k*T/elc)**2*32*pi*a*tanh(psi_s_absence/4)**2*exp(-kappa*h_a)
    V_vdW_max_a = -A*a/(12*h_a)
    Vmax_a = V_EL_max_a + V_vdW_max_a
    print("V max (absence) = ", Vmax_a/(k*T))

    # real aggregation rate constant 
    if W_50>1:
        kstar_50 = k_peri/W_50 
    else:
        kstar_50 = k_peri

    # aggregation time 
    print("tau (50)    =", 2/(kstar_50*n), 's' )
    
    ## Presence of pregnancy
    print("\n--------------------------------------------")
    print("I (M = mol/L)=", I_list[i], "and Vmax_25 = ", Vmax_presence[i])

    # Calculate the stability ratio
    W_25 = (1/(kappa*2*a))*exp(Vmax_presence[i])
    print("W_25 (presence) = ", W_25)
    W_25_list.append(W_25)
    
    # Ionic strength
    F_EL_0_p = eps * (k*T/elc)**2 * 32 * pi *kappa* a * tanh(psi_s_presence/4)**2
    F_vdw_0_p = - A * a /12
    res = fsolve( lambda h_p: F_EL_0_p*exp(-kappa*h_p) + F_vdw_0_p/h_p**2, 2.0e-9 )
    # fsolve returns an array with 1 value. The following line is just to extract this value
    h_p = res[0]
    print("h =", h_p )
    
     # Define the potentials at the location of the max h 
    V_EL_max_p = eps*(k*T/elc)**2*32*pi*a*tanh(psi_s_presence/4)**2*exp(-kappa*h_p)
    V_vdW_max_p = -A*a/(12*h_p)
    Vmax_p = V_EL_max_p + V_vdW_max_p
    print("V max (presence)= ", Vmax_p/(k*T))
    
    # CCC presence
    gamma_0_presence = tanh(psi_s_presence/4.)
    print("gamma_0_presence = ", gamma_0_presence)
    CCC_presence = 3.8e-36*(gamma_0_presence**4)/(A**2*z**6)
    print("CCC presence (mol/m3) = ", CCC_presence)

    # real aggregation rate constant 
    if W_25>1:
        kstar_25 = k_peri/W_25
    else:
        kstar_25 = k_peri

    # aggregation time 
    print("tau (25)    =", 2/(kstar_25*n), 's' )


    print("W_50_list = ", W_50_list)
    print("W_20_list = ", W_25_list)

# Plot the stability ratio
plt.figure(0)
plt.scatter(I_list, W_50_list, label = "Vmax_50 (absence)", color = "red", marker = "o")
plt.scatter(I_list, W_25_list, label = "Vmax_20 (presence)", color = "blue", marker = "o")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("I (M = mol/L)")
plt.ylabel("W")
plt.title("Stability ratio")
plt.legend()
plt.grid()
plt.savefig("pregnancy_w.png")
