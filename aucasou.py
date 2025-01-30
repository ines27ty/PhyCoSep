from math import *
from pylab import *
from ct import *            # File of constants
import numpy as np  
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


a = 25e-9          # rayon des particules (m)
psi_s = 0.040 * elc/kT            # dimensionless potential  psi = Psi * e / kT
phi = 0.00232         # fraction volumique en particules
A = 2e-21           # constante Hamaker silie - eau (J)
T = 293
eta   = 0.001                     # Pa.s

## Partie 2
print("\n----------PARTIE 2---------------------------------")
print("Nouvelle valeur de phi : ", phi)


# Q1 : Surface totale des particules
print("\n--------------------------------------------")
print("Q1 (partie 2)")

V = 1       # volume m3
n = 3*V*phi/(4*np.pi*a**3)
print("n = ", n)
S = n*4*np.pi*a**2

print("S (m2)=", S)

# Q2 : 
print("\n--------------------------------------------")
print("Q2 (partie 2)")
print("Nombre de charges de surface de 1m3 : n (nombre/m3) = ", n)
n0 = phi / (4/3*np.pi*a**3)
print("Nombre de particules par m3 : n0 (nombre/m3) = ", n0)

sigma_p = 0.5       # e/nm
c_s = sigma_p*S*1e18*0.001/Na
print("cs (mol/L) = ", c_s)

# Q3 : 
print("\n--------------------------------------------")
print("Q3 (partie 2)")
lambda_d = np.sqrt(eps*kT/(2*1000*Na*c_s*elc**2))
print("lambda_d (m) = ", lambda_d)
kappa = 1/lambda_d
print("kappa (m-1) =", kappa)
ka = kappa *a 
print("ka=", ka, ">1 donc Derjaguin valable mais limite, donc approximation LSA valable aussi sûrement")


## Q4 : Derjaguin approximation for the electrostatic and vdW forces, valid if ka>>1 and h<<a 
print("\n----------Derjaguin POUR TESTER si valide ou pas ----------------------------------")
print("Q4")
r = np.linspace((2*a + 0.01e-9), (2*a + 15e-9), 100)
h= r - 2*a
V_EL = eps * (kT/elc)**2*32*pi*a*tanh(psi_s/4)**2*exp(-kappa*(h))/(kT)
V_vdw = - A*a/(12*h)/(kT)
V_DLVO = V_EL + V_vdw           # Potentiel d'intéraction total DLVO par unité de kT

plt.figure(0)
plt.plot(h,V_EL, label="V_EL")
plt.plot(h,V_vdw, label="V_vdw")
plt.plot(h,V_DLVO, label="V_DLVO")
plt.xlabel("h (m)")
plt.ylabel("V")
plt.ylim(-15,40)
plt.legend()
plt.grid()
plt.savefig("potentiel_partie2_Derjaguin.png")

## Q6 : 
print("\n--------------------------------------------")
print("Q6")
V_DLVO_max = max(V_DLVO)
print("V_DLVO_max/(kT) = ", V_DLVO_max, ">> 15 donc stable. On peut garder le même lot de Klebosol")
# Comparer à kT

##Q8 : Rapport de la stabilité
print("\n--------------------------------------------")
print("Q8")
W = 1/(2*a*kappa)*exp(V_DLVO_max)
print("W = ", W)

print("\n---------------STOP -----------------------------")


##Q9 : 
print("\n--------------------------------------------")
print("Q9")
k_peri = 8*kT/(3*eta)
print("Taux d'aggrégation (sans correction) en supposant péricinétique : k_peri = ", k_peri)
if W>1:
    kstar = k_peri/W 
else:
    kstar = k_peri
print("Taux d'aggrégation (avec correction) : kstar = ", kstar)

print("Temps d'aggrégation : tau = ", 2/(kstar*n), "s")
#en minutes, heures, jours
print("Temps d'aggrégation : tau = ", 2/(kstar*n)/60, "min")
print("Temps d'aggrégation : tau = ", 2/(kstar*n)/3600, "h")
print("Temps d'aggrégation : tau = ", 2/(kstar*n)/3600/24, "jours")

# Q10 
print("\n--------------------------------------------")
print("Q10")
print("Comme on a un grand temps d'aggrégation : tau = ", 2/(kstar*n)/3600/24, "jours, alors on peut dire que le système est stable en conditions de stockage.")

