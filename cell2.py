#!/usr/bin/env python
import os, sys
from pylab import *
from scipy.optimize import fsolve, minimize
from scipy.integrate import simps
#sys.path.append(os.path.abspath('/home/hallez/pobos_svn/tools/prep/'))
#from prep import *
import ct
import mpmath
ion()


class CellModelSimulation():
  def __init__(self):
    self.adim  = False 
    self.eps   = 78.
    self.a     = -1. #nm
    self.I     = -1.      # mol/L
    self.sigma = -1.    # e/nm**2
    self.phi   = -1. 
    self.BC    = "FixedCharge"  # FixedCharge or FixedPotential
    self.lin   = False
    self.quiet = True
    self.saltfree = False
    self.resolution = 400
    self.n_0   = -1.
    self.Pi    = 0.
    self.gamma0 = 0.
    self.psi    = 0.


  def SetParameters(self,eps=78, a=-1, I=-1, sigma=-1, phi=-1, BC="FixedCharge", pH=-1, pK=[], Gamma_tot=[], CS=[], ksi=[], lin=False, quiet=True, saltfree=False, res=400, adim=False):
    self.eps = eps
    self.a = a
    self.I = I
    self.sigma = sigma
    self.phi = phi
    self.BC = BC
    self.pH = pH
    self.pK = pK
    self.Gamma_tot = Gamma_tot
    self.CS = CS
    self.ksi = ksi
    self.rho = self.phi / (4./3.*pi*(self.a*1.0e-9)**3)
    self.lin = lin
    self.quiet = quiet
    self.saltfree = saltfree
    self.resolution = res
    self.adim = adim
    self.ComputeConstants()
    self.BuildMesh()


  def ComputeConstants(self):
    # Check 
    if self.eps<0 or self.a<0 or self.I<0 or self.phi<0:
      print("Parameters need being set: eps, a, I, sigma, phi")
      sys.exit()
    # Bulk number density:
    self.n_0 = self.I*1000*ct.Na
    self.n_c = 0.
    S = 4*pi*self.a**2  # surface in nm**2
    # Compute initial surface charge from the planar solution 
    if not self.pH==-1:
      res        = fsolve( lambda psiD: self.compute_sigmaCR(psiD)-2.0*sinh(psiD/2.), 0.0, xtol=1.0e-12)
      psiD       = res[0]
      self.sigma = 2.0*sinh(psiD/2.0)
    self.Z = self.sigma * S  # nb charges
    if self.saltfree:
      # in this case the actual nb density scale is not 2n_salt but n_counterion=Z*phi/(1-phi)/Vcoll
      self.n_c = self.Z*self.phi/(1-self.phi)/(4/3.*pi*(self.a*1.0e-9)**3)
      self.n_0 = 0.
    # Debye length
    self.lambdaD = sqrt(ct.eps0*self.eps*ct.kT/((2*self.n_0+self.n_c)*ct.elc**2))
    # Non-dimensional particle radius (a comes in nm)
    self.a = self.a*1e-9/self.lambdaD
    # Radius of the cell (non-dimensional)
    self.R = self.a/self.phi**(1./3.)
    self.V = 4./3.*pi*self.R**3
    # Non-dimensional charge (sigma comes in e/nm**2)
    self.sigmascale = sqrt( (2*self.n_0+self.n_c)*ct.eps0*self.eps*ct.kT)
    if self.pH==-1:
      self.sigma     = self.sigma    *ct.elc*1.0e18/self.sigmascale
    else:
      self.Gamma_tot = self.Gamma_tot*ct.elc*1.0e18/self.sigmascale
      # Stern capacity
      self.CSscale = self.sigmascale * ct.elc/ct.kT
      self.CS      = self.CS / self.CSscale
    # Bjerrum length
    self.lambdaB = ct.elc**2/(4*pi*ct.kT*ct.eps0*self.eps)

  def PrintConstants(self):
    print("Debye length : ", self.lambdaD)
    print("Radius : ", self.a)
    print("Surface charge density : ", self.sigma)
    print("Volume fraction : ", self.phi)
    print("Cell radius : ", self.R)

  def BuildMesh(self):  
    self.dr = (self.R-self.a)/self.resolution;
    self.r = arange(self.a, self.R, self.dr)
    return

  def f_const_charge(self, psi):
    # Update charge density if necessary
    if not self.pH==-1:
      self.compute_sigmaCR(psi[0])
    out = zeros(len(self.r))
    if not self.lin:
      out[0]    = sinh(psi[   0]) - \
	          1./self.r[0   ]**2*( (self.r[   0]+self.dr*0.5)**2*(psi[1   ]-psi[0   ])/self.dr - \
		                       (self.r[   0]-self.dr*0.5)**2*(psi[0   ]-psi[1   ]-self.sigma*2*self.dr)/self.dr \
				     )/self.dr
      out[1:-1] = sinh(psi[1:-1]) - \
	          1./self.r[1:-1]**2*( (self.r[1:-1]+self.dr*0.5)**2*(psi[2:  ]-psi[1:-1])/self.dr - \
		                       (self.r[1:-1]-self.dr*0.5)**2*(psi[1:-1]-psi[0:-2])/self.dr \
				     )/self.dr
      out[-1]   = sinh(psi[  -1]) - \
	          1./self.r[  -1]**2*( (self.r[  -1]+self.dr*0.5)**2*(psi[  -2]-psi[  -1])/self.dr - \
		                       (self.r[  -1]-self.dr*0.5)**2*(psi[  -1]-psi[  -2])/self.dr \
				     )/self.dr
      if self.saltfree:
        out = zeros(len(self.r))
        # NON DIM CORRECTLY HERE : REGARDE DANS TON CAHIER !
        out[   0] = exp(psi[   0]) - \
	            1./self.r[   0]**2*( (self.r[   0]+self.dr*0.5)**2*(psi[1   ]-psi[0   ])/self.dr - \
	                                 (self.r[   0]-self.dr*0.5)**2*(psi[0   ]-psi[1   ]-self.sigma*2*self.dr)/self.dr \
	                               )/self.dr
        out[1:-1] = exp(psi[1:-1]) - \
	            1./self.r[1:-1]**2*( (self.r[1:-1]+self.dr*0.5)**2*(psi[2:  ]-psi[1:-1])/self.dr - \
		                         (self.r[1:-1]-self.dr*0.5)**2*(psi[1:-1]-psi[0:-2])/self.dr \
				       )/self.dr
        out[  -1] = exp(psi[  -1]) - \
	            1./self.r[  -1]**2*( (self.r[  -1]+self.dr*0.5)**2*(psi[  -2]-psi[  -1])/self.dr - \
		                         (self.r[  -1]-self.dr*0.5)**2*(psi[  -1]-psi[  -2])/self.dr \
				       )/self.dr
    else:

      out = zeros(len(self.r))
      out[   0] = self.gamma0 + psi[   0] - \
	          1./self.r[   0]**2*( (self.r[   0]+self.dr*0.5)**2*(psi[   1]-psi[   0])/self.dr - \
		                       (self.r[   0]-self.dr*0.5)**2*(psi[   0]-psi[   1]-self.sigma*2*self.dr)/self.dr \
				     )/self.dr
      out[1:-1] = self.gamma0 + psi[1:-1] - \
	          1./self.r[1:-1]**2*( (self.r[1:-1]+self.dr*0.5)**2*(psi[2:  ]-psi[1:-1])/self.dr - \
		                       (self.r[1:-1]-self.dr*0.5)**2*(psi[1:-1]-psi[0:-2])/self.dr \
				     )/self.dr
      out[  -1] = self.gamma0 + psi[  -1] - \
      	          1./self.r[  -1]**2*( (self.r[  -1]+self.dr*0.5)**2*(psi[  -2]-psi[  -1])/self.dr - \
      		                       (self.r[  -1]-self.dr*0.5)**2*(psi[  -1]-psi[  -2])/self.dr \
      				     )/self.dr
      #out[  -1] = - (4*psi[-2]-psi[-3]-3*psi[-1]) #/(2*self.dr)

    return out

  def f_const_pot(self,psi):
    out = zeros(len(self.r))
    out[   0] = psi[0]-self.pot
    out[1:-1] = sinh(psi[1:-1]) - \
	        1./self.r[1:-1]**2*( (self.r[1:-1]+self.dr*0.5)**2*(psi[2:  ]-psi[1:-1])/self.dr - \
		                     (self.r[1:-1]-self.dr*0.5)**2*(psi[1:-1]-psi[0:-2])/self.dr \
				   )/self.dr
    out[  -1] = sinh(psi[  -1]) - \
	        1./self.r[  -1]**2*( (self.r[  -1]+self.dr*0.5)**2*(psi[  -2]-psi[  -1])/self.dr - \
		                     (self.r[  -1]-self.dr*0.5)**2*(psi[  -1]-psi[  -2])/self.dr \
				   )/self.dr
    return out


  def solve(self):
    
    if self.BC=="FixedCharge":
      # Estimate of the surface potential from the surface charge density
      #psi_s, charge_verif = charge_to_pot(self.sigma, self.a, -100., 100.);
      #if not self.quiet:
      #  print "psi_s, charge_verif, sigma", psi_s, charge_verif, sigma
      
      #if abs(self.sigma-charge_verif)>1.0e-10:
      #  print "error in computing surface potential"
      #  sys.exit(0)
      # Launch computations with constant charge condition
      initpsi = 0.
      ier = 0
      while not ier == 1:
        initpsi += 1
        print(initpsi)
        psi_0 = initpsi*ones(len(self.r))
        #psi_0    = ones(len(self.r))*psi_s
        x, infodict, ier, mesg = fsolve( lambda pot: self.f_const_charge(pot), psi_0,full_output=True,xtol=1.0e-12, maxfev=100000000)
      self.psi = x
      if not ier==1:
          print("solution did not converge")
          print(mesg)
          sys.exit()
    elif self.BC=="FixedPotential":
      # Launch computations with constant potential condition
      psi_0    = ones(len(self.r))*psi_s/2.
      self.psi = fsolve( lambda pot: f_const_pot(self,pot),psi_0,xtol=1.0e-12, maxfev=100000000)
    else:
      print("Unknown boundary condition. Available ones are FixedCharge or FixedPotential\n")
      sys.exit()


  def ExactSolutionDHCM(self):
    rap = (self.R-1)/(self.R+1)
    # potential
    alpha = self.sigma*self.a**2 / ((1-self.a)*exp(self.a)+(1+self.a)*rap*exp(2*self.R-self.a))
    beta  = rap*exp(2*self.R)
    self.psi = alpha /self.r * (exp(self.r)+beta*exp(-self.r)) \
    # osmotic pressure	       
    self.Pi = self.psi[-1]**2/2.     
    if not self.adim:
      if self.saltfree:
        self.Pi = self.Pi * self.n_c*ct.kT
      else:
        self.Pi = self.Pi * 2*self.n_0*ct.kT
    # free energy
    self.A = -2*pi*alpha**2*(exp(2*self.R)-exp(2*self.a)-beta**2*(exp(-2*self.R)-exp(-2*self.a))) + float( (4*pi*alpha**2*( beta*(1./self.R-1./self.a) -mpmath.gammainc(0,-2*self.R,-2*self.a) - mpmath.gammainc(-1,-2*self.R,-2*self.a) \
	                -beta**2*(   mpmath.gammainc(0,2*self.a,2*self.R) + mpmath.gammainc(-1,2*self.a,2*self.R) ) ) + 4.0*pi*self.a**2*self.sigma*self.psi[0]).real) 
    if not imag(self.A)==0:
      print("Stopping ! free energy is complex in ExactSolutionDHCM() !")
      sys.exit(0)
    self.A = real(self.A)
    # dimensions
    self.A *= 2*self.n_0*ct.kT*self.lambdaD**3

  def ExactSolutionDHCM_lin_av_pot(self):
    psi_av = -arcsinh(3*self.Z/(8*pi*self.n_0*(self.R**3-self.a**3)*self.lambdaD**3))
    kappa  = 1./self.lambdaD*sqrt(cosh(psi_av))
    n_plus  = self.n_0*exp(-psi_av)
    n_minus = self.n_0*exp(+psi_av)
    a       = self.a*self.lambdaD
    r       = self.r*self.lambdaD 
    R       = self.R*self.lambdaD
    rap     = (kappa*R-1)/(kappa*R+1)
    A1 = self.Z*ct.lambdaB / ((kappa*a-1.0)*exp(kappa*a)-(1+kappa*a)*rap*exp(2*kappa*R-kappa*a))
    A2 = A1*rap*exp(2*kappa*R)
    A3 = psi_av + (n_plus-n_minus)/(n_plus+n_minus)
    self.psi = A1*exp(kappa*r)/r + A2*exp(-kappa*r)/r+A3
    # LCM osmotic pressure according to Deserno
    #self.Pi = self.psi[-1]**2/2.     
    self.Pi = cosh(psi_av) + sinh(psi_av)*(self.psi[-1]-psi_av) + 0.5*cosh(psi_av)*(self.psi[-1]-psi_av)**2  -1.0   
    if not self.adim:
      self.Pi = self.Pi * 2*self.n_0*ct.kT
    # linearized microion grand potential per macroion
    Vprime = 4/3.*pi*(R**3-a**3)
    N_plus  = n_plus  * Vprime
    N_minus = n_minus * Vprime
    x_plus  = N_plus  / 1.0
    x_minus = N_minus / 1.0
    omega_lin = x_plus*(log(n_plus/self.n_0)-1.0) + x_minus*(log(n_minus/self.n_0)-1.0) + 0.5*self.Z*(psi_av-self.psi[0])
    print(omega_lin,"aa")
    ## free energy
    #self.A = -2*pi*alpha**2*(exp(2*self.R)-exp(2*self.a)-beta**2*(exp(-2*self.R)-exp(-2*self.a))) + float( (4*pi*alpha**2*( beta*(1./self.R-1./self.a) -mpmath.gammainc(0,-2*self.R,-2*self.a) - mpmath.gammainc(-1,-2*self.R,-2*self.a) \
    #                    -beta**2*(   mpmath.gammainc(0,2*self.a,2*self.R) + mpmath.gammainc(-1,2*self.a,2*self.R) ) ) + 4.0*pi*self.a**2*self.sigma*self.psi[0]).real) 
    #if not imag(self.A)==0:
    #  print "Stopping ! free energy is complex in ExactSolutionDHCM() !"
    #  sys.exit(0)
    #self.A = real(self.A)
    ## dimensions
    #self.A *= 2*self.n_0*ct.kT*self.lambdaD**3

  def ComputePressure(self):
    # compute dimensionless osmotic pressure
    if (not self.lin) and (not self.saltfree):
      self.Pi =  cosh(self.psi[-1]) - 1.
    elif self.lin:
      self.Pi =  (self.psi[-1]+arctanh(self.gamma0))**2/2.
    elif self.saltfree:
      self.Pi = self.n_c*exp(self.psi[-1])
    # get dimensions back
    if not self.adim:
      if self.saltfree:
        self.Pi = self.Pi * self.n_c*ct.kT
      else:
        self.Pi = self.Pi * 2*self.n_0*ct.kT

  def CountIons(self): 
    if (not self.lin) and (not self.saltfree):
      Nm = trapz( self.n_0*exp(-self.psi)*4*pi*(self.r*self.lambdaD)**2, self.r*self.lambdaD)
      Np = trapz( self.n_0*exp(+self.psi)*4*pi*(self.r*self.lambdaD)**2, self.r*self.lambdaD)
      N_col = self.Z
    return N_col, Np-N_col, Nm

  
  def ComputeVirialPressure(self):
    if (not self.lin) and (not self.saltfree):
      # (Eself is not analytical for a sphere IN THE NON LINEAR REGIME)X 
      # -> WELL ! it seems it is actually the same ! see PoBoS doc !
      # but using the analytical value prevents relying on numerical error cancellation. 
      Eself = 2*pi*self.a**3*self.sigma**2 # 0.0 #314.073  
      P_Vion = trapz( (cosh(self.psi)-1)*4*pi*self.r**2 ,x=self.r ) / self.V
      P_Sion = 4*pi*self.a**2/(3*self.V) * (cosh(self.psi[0])-1) * self.a 
      P_ESion = -1./(3*self.V)*trapz( self.psi*sinh(self.psi)/2. *4*pi*self.r**2 ,x=self.r )
      P_EScol =  4*pi*self.a**2/(3*self.V) * 0.5*self.sigma*self.psi[0] 
      P_self  = Eself / (3*self.V)
      self.Pi = P_Vion + P_Sion + P_ESion + P_EScol - P_self 

      self.ComputeElectricField()
      print(P_Vion, P_Sion, P_ESion, P_EScol, -P_self)
      print("debug", trapz( cosh(self.psi)*4*pi*self.r**2 ,x=self.r ) / self.V)
      #print "debug", cosh(self.psi[-1])-(self.a/self.R)**3*cosh(self.psi[0])
      #print "debug", -trapz(self.r*(-self.E)*sinh(self.psi)*4*pi*self.r**2,x=self.r)/(3*self.V)
      #print "debug", cosh(self.psi[-1])-(self.a/self.R)**3*cosh(self.psi[0])-trapz(self.r*(-self.E)*sinh(self.psi)*4*pi*self.r**2,x=self.r)/(3*self.V)
      print("debug2", trapz(self.r*(-self.E)*sinh(self.psi)*4*pi*self.r**2,x=self.r)) #/ (3*self.V)
      print("debug2", trapz( 0.5*self.E**2 *4*pi*self.r**2 ,x=self.r )-2*pi*self.a**3*self.sigma**2)
    elif self.lin:
      # lin for now
      # Eself is analytical for a sphere IN THE LINEAR REGIME 
      Eself = 2*pi*self.a**3*self.sigma**2
      P_Vion = trapz( (1+self.psi**2/2.-1)*4*pi*self.r**2 ,x=self.r ) / self.V
      P_Sion = 4*pi*self.a**2/(3*self.V) * (1+self.psi[0]**2/2.-1) * self.a 
      P_ESion = -1./(3*self.V)*trapz( self.psi**2/2. *4*pi*self.r**2 ,x=self.r )
      P_EScol =  4*pi*self.a**2/(3*self.V) * 0.5*self.sigma*self.psi[0] 
      P_self  = Eself / (3*self.V)
      self.Pi = P_Vion + P_Sion + P_ESion + P_EScol - P_self 
    elif self.saltfree:
      pass
    # get dimensions back
    if not self.adim:
      if self.saltfree:
        pass
      else:
        self.Pi = self.Pi * 2*self.n_0*ct.kT

  def ComputeElectricField(self):
    self.E = zeros(len(self.r))
    dr = self.r[1]-self.r[0]
    self.E[0]    = - (  (4*self.psi[1]-self.psi[2]-3*self.psi[0])/(2.0*dr)       )
    self.E[1:-1] = - (  (self.psi[2:]-self.psi[:-2])/(2*dr)                      )
    self.E[-1]   = - ( -(4*self.psi[-2]-self.psi[-3]-3*self.psi[-1])/(2.0*dr)    )

  def ComputeFreeEnergy(self):
    # volume integral
    self.ComputeElectricField()
    if self.lin:
      self.A = simps((-self.psi**2/2.-0.5*self.E**2)*(4*pi*self.r**2),x=self.r)
    else:
      self.A = simps((ones(len(self.r))-cosh(self.psi)-0.5*self.E**2)*(4*pi*self.r**2),x=self.r)
    # surface integral
    self.A += 4*pi*self.a**2*self.sigma*self.psi[0]
    # dimensions
    self.A *= 2*self.n_0*ct.kT*self.lambdaD**3

  def ComputeGrandPotential(self):
    # volume integral
    n_plus  = exp(-self.psi)  #dimensionless
    n_minus = exp(+self.psi)
    if self.lin:
      print('error')
      sys.exit(0)
    else:
      self.Omega = trapz( (-cosh(self.psi)+self.psi*sinh(self.psi)*0.5 +1.0)*(4*pi*self.r**2),x=self.r)
      # +1.0 is a trick!!!
    # surface integral
    self.Omega += 0.5 * 4*pi*self.a**2*self.sigma*self.psi[0]
    # dimensions
    self.Omega *= 2*self.n_0*ct.kT*self.lambdaD**3 

  def GrandPotential_PBDFT(self,psi,r,a,sigma):
    n_plus  = exp(-psi)  #dimensionless
    n_minus = exp(+psi)
    Omega = trapz( (-cosh(psi)+psi*sinh(psi)*0.5)*4*pi*r**2,x=r)
    # surface integral
    Omega += 0.5 * 4*pi*a**2*sigma*psi[0]
    return abs(Omega)


  def test(self):
    dr = self.r[1]-self.r[0]
    #dpsidr = zeros(len(self.psi))
    #dpsidr[0] = (self.psi[1]-self.psi[0])/dr
    #dpsidr[1:-1] = (self.psi[2:]-self.psi[:-2])/(2*dr)
    #dpsidr[-1] = (self.psi[-1]-self.psi[-2])/dr
    self.ComputeElectricField()
    V = 4./3.*pi*(self.r[-1]**3-self.r[0]**3) #adim 
    # a la belloni, zi.e.E = force sur 1 ion, c est une pression (pas osmo)
    integ = simps( self.E*sinh(self.psi)*self.r**3*4*pi, x=self.r) / (3*V) 
    print("test", integ*2*self.n_0*ct.kT, 2*self.n_0*ct.kT, self.r[0]*ct.kT/self.lambdaD**3*self.E[0]*self.Z/(3*V))

  def renormalize(self,type='SCR'):
    """
    Renormalize following Trizac et al. (2003)'s Alexander's prescription revisited.
    """
    # Solve the non-linear problem for constant charge 
    self.lin = False
    self.solve()
    # Get the potential at the cell boundary
    psi_R = self.psi[-1]
    # Compute the inverse screening length at cell boundary 
    self.kappa_PB = sqrt(cosh(psi_R))/self.lambdaD
    self.gamma0 = tanh(psi_R)
    if type=='SCR':
      # Compute effective charge from eq (16) of Trizac et al. 2003
      self.Zeff = self.gamma0/(self.kappa_PB*self.lambdaB)*( (self.kappa_PB**2*self.a*self.R*self.lambdaD**2-1.0)*sinh(self.kappa_PB*(self.R-self.a)*self.lambdaD)+ \
                                                    self.kappa_PB*(self.R-self.a)*self.lambdaD*cosh(self.kappa_PB*(self.R-self.a)*self.lambdaD) ) 
    elif type=='EPC':
      kR = self.kappa_PB*self.R*self.lambdaD
      ka = self.kappa_PB*self.a*self.lambdaD
      self.Q = self.gamma0/(self.kappa_PB*self.lambdaB)*( kR*cosh(kR)-sinh(kR))
      self.Zeff = (1+ka)*exp(-ka)*self.Q
      print(self.Q, self.Zeff)
    else: 
      print("Unknown renormalization type. Try nothing or EPC")
    # Compute effective salt number density (useless if just renormalizing)
    V  = 4./3.*pi*(self.R*self.lambdaD)**3
    nc = 1.0 / V
    self.nseff = 0.5*( self.kappa_PB**2/(4*pi*self.lambdaB)*(1.0-self.gamma0**2)*(1.0-self.phi) - self.Zeff*nc*(1.0-self.gamma0) )  # a priori useless
    self.n0eff = ct.eps0*self.eps*ct.kT/(2.*ct.elc**2*(1./self.kappa_PB)**2)
    self.sigma_eff_adim = self.Zeff*ct.elc/(4*pi*(self.a*self.lambdaD)**2) / sqrt( 2*self.n0eff*ct.kT*ct.eps0*self.eps )
    print("Original  parameters:")
    print("phi            =", self.phi)
    print("Z              =", self.Z)
    print("n0             =", self.n_0)       
    print("a    adim      =", self.a)
    print("sigma     adim =", self.sigma)
    print("Effective parameters:")
    print("Zeff           =", self.Zeff)
    print("kappa_PB       =", self.kappa_PB)
    print("gamma0         =", self.gamma0)
    print("n0eff          =", self.n0eff)
    print("nseff          =", self.nseff)
    print("aeff adim      =", self.a*self.lambdaD*self.kappa_PB)
    print("sigma_eff adim =", self.sigma_eff_adim) 

  def compute_sigmaCR(self,psiD):
    # compute a charge value based on chemistry and the surface potential psiD
    x, infodict, ier, mesg = fsolve(lambda sigmaI: sigmaI - sum( self.Gamma_tot *( 1./(1.+10.0**(self.pH-self.pK)*exp(psiD + sigmaI/self.CS)) + self.ksi) ), 0., xtol=1.0e-12, full_output=True)
    self.sigma = x[0]
    if not ier==1:
      print("solution did not converge in compute_sigmaCR")
      print(mesg)
      sys.exit()
    return self.sigma

  def TestFreeEnergy(self):
    phi_list = arange(0.02,0.221, 0.02)
    A_num    = []
    A_theo   = []
    close('all')
    figure()
    self = CellModelSimulation()
    for phi in phi_list:
      self.SetParameters(a=7.5,I=.001559,res=400,phi=phi,sigma=0.0366,lin=True)
      self.ExactSolutionDHCM()
      A_theo.append(self.A)
      self.ComputeFreeEnergy()
      A_num.append(self.A)
    plot(phi_list, array(A_num)/ct.kT, 'o')
    plot(phi_list, array(A_theo)/ct.kT, 'r')
    xlabel('$\phi$')
    ylabel('Free Energy (kT)')
    legend(('Numerical integration','Theoretical solution'), loc='lower right')
    
  def Plot2D(self):
    fig, ax = subplots(subplot_kw=dict(projection='polar'))
    ntheta = 1000
    theta   = linspace(0,360,ntheta)
    psi_for_plot = zeros( (len(self.r)+1,ntheta) )
    for itheta in range(ntheta):
      psi_for_plot[1:,itheta] = self.psi
    psi_for_plot[0,:] = self.psi[0]*1000000000.0  
    CF = ax.contourf(theta,concatenate( ([0],self.r) ),psi_for_plot,linspace(min(self.psi),max(self.psi),30))
    CF.set_clim(0.,1.5934) #min(self.psi),max(self.psi))
    CF.cmap.set_over('w')
    colorbar(CF)
    axis('off')
    savefig('cell_a2.png')
    #from matplotlib.patches import Polygon
    #from matplotlib.collections import PatchCollection
    #circtab = zeros((ntheta,2))
    #circtab[:,0] = linspace(0,360,ntheta)
    #circtab[:,0] = self.a
    #circle = Polygon( circtab, True )
    #patches = [circle]
    #p = PatchCollection(patches,color='white')
    #ax.add_collection(p)



  
  def TestLin(self):
    close('all')
    figure()
    self = CellModelSimulation()
    #self.SetParameters(a=5,I=0.02259,res=400,phi=0.20,sigma=0.05) # 0.0366
    #self.SetParameters(a=50.0,I=.0001,res=400,phi=0.20,sigma=0.015915494309189534,lin=True)
    self.SetParameters(a=7.5,I=.001558*4.4,res=400,phi=0.15,sigma=0.5/2.67,lin=True)
    self.PrintConstants()
    self.solve()
    self.ComputePressure()
    print("Lin cell model pressure (numerical)",self.Pi)
    plot(self.r, self.psi, 'o')
    self.ExactSolutionDHCM()
    print("Lin cell model pressure (analytical)",self.Pi)
    plot(self.r, self.psi, '-g')
    self.ComputeVirialPressure()
    print("Lin cell model pressure (virial/analytical sol)",self.Pi)
    if self.adim:
      print("Lin cell model Eself                           ",self.Pi * 3 * self.V) 
    else:
      print("Lin cell model Eself                           ",self.Pi * 3 * self.V * self.lambdaD**3) 

    self.ExactSolutionDHCM_lin_av_pot()
    print("Lin cell model pressure (analytical, lin around av pot)",self.Pi)
    plot(self.r, -self.psi, '-r')
    self.ComputePressure()
    print("Lin cell model Helmholtz free energy ", self.A/ct.kT)
    self.ComputeFreeEnergy()
    print("Numerically obtained Helmholtz free energy ", self.A/ct.kT)
    #compute pressure from 2 values of free energy
    self.SetParameters(a=7.5,I=.001559,res=400,phi=0.1601,sigma=0.0366,lin=True)
    self.ExactSolutionDHCM()
    A2 = self.A
    self.SetParameters(a=7.5,I=.001559,res=400,phi=0.1599,sigma=0.0366,lin=True)
    self.ExactSolutionDHCM()
    A1 = self.A
    Pi = (A2-A1)/(.1601-.1599)*4./3*pi*self.a**3/(4./3*pi*self.R**3)**2 / self.lambdaD**3
    print(self.a)
    print("Cell model pressure by differentiating theo free energy ", Pi)



  def TestNL(self):
    close('all')
    #figure()
    self = CellModelSimulation()
    #self.SetParameters(a=5,I=0.02259,res=400,phi=0.20,sigma=0.05) # 0.0366
    self.SetParameters(a=7.5,I=0.001558,res=800,phi=0.1,sigma=0.5,adim=True) # 0.0366
    #self.SetParameters(a=50.,I=0.0002,res=400,phi=0.08,sigma=0.015915494309189534)
    #self.SetParameters(a=73.0,I=4.239e-6,res=400,phi=0.01,sigma=1.5)
    self.PrintConstants()
    # compute linear solution
    self.ExactSolutionDHCM()
    plot(self.r/self.a, self.psi, '--b', lw=2)
    # compute non linear solution
    self.solve()
    plot(self.r/self.a, self.psi, '-b', lw=2)
    self.ComputePressure()
    print("NL cell model pressure ",self.Pi)
    self.ComputeVirialPressure()
    print("Virial Pressure: ", self.Pi)
    print("Eself", 2*pi*self.a**3*self.sigma**2)
    print("Eself vir", self.Pi*3*4/3.*pi*self.R**3)
    self.ComputeGrandPotential()
    self.test()
    print("NL Grand Potential : ", self.Omega)
    print("NL pot at R : ", self.psi[-1])
    self.renormalize()
    psiNL = self.psi

    #compute net charge as a function of r
    import scipy.integrate as integ
    Q = integ.simps(2*self.n_0*cosh(self.psi[:161])*4*pi*(self.r[:161]*self.lambdaD)**2,x=self.r[:161]*self.lambdaD)
    print("Q=", Q)

 
    # test DFT
    cons = ({'type':'eq','fun':lambda psi: self.Z - 2*self.n_0*self.lambdaD**3*trapz(sinh(psi)*4*pi*self.r**2,self.r)},\
            {'type':'eq','fun':lambda psi: (psi[1]-psi[0])/self.dr+self.sigma})
    res = minimize(lambda psi: self.GrandPotential_PBDFT(psi,self.r,self.a,self.sigma), 3.*ones(len(self.r)), method='SLSQP',constraints=cons,tol=1e-8)
    print(res)
    self.psi = res.x
    plot(self.r/self.a,self.psi,'-+')




    Ieff = self.n0eff/(1000*ct.Na)
    sigmaeff = self.Zeff*ct.elc/(4*pi*(self.a*self.lambdaD)**2) / ct.elc/1.0e18
    self.gamma0 = tanh(arccosh(self.Pi/(2*self.n_0*ct.kT)+1)) #0.974305186225 a 15%
    print(Ieff, sigmaeff)
    self.SetParameters(a=7.5,I=Ieff,res=400,phi=0.01,sigma=sigmaeff,lin=True ) # 0.0366
    self.PrintConstants()
    self.solve()
    plot(self.r/self.a, self.psi+arctanh(self.gamma0), '-r', lw=2)

    ylim(0,10)
    xlim(1,4.5)
    rcParams['text.usetex']=True
    rcParams['font.size']=16   
    legend( ('LPB','PB',"LPB renormalis\\'e"), loc='upper right')
    xlabel('$r/a$',fontsize=22)
    ylabel('$\\tilde{\\psi}$',fontsize=22)

    ###find pos of Rep
    ##sigmael = 0.056
    ##phi = 0.01
    ##R = (self.a*self.lambdaD)*phi**(-1./3)
    ##fp = (self.kappa_PB*R+1)/(2*self.kappa_PB)*exp(-self.kappa_PB*R)
    ##fm = (self.kappa_PB*R-1)/(2*self.kappa_PB)*exp(+self.kappa_PB*R)
    ##res = fsolve(lambda Rep:self.gamma0*(fp*exp(self.kappa_PB*Rep)*(Rep*self.kappa_PB-1)/Rep**2+fm*exp(-self.kappa_PB*Rep)*(-Rep*self.kappa_PB-1)/Rep**2)+sigmael*ct.elc*1.0e18/ct.eps * ct.elc/ct.kT,self.a*self.lambdaD*2)
    ##print "res (rayon electroph)", res
    ###print [(i,self.r[i]*self.lambdaD) for i in range(len(self.r))]
    ##plot([res, res],[0,max(self.psi)],'k')
    ###print "psi'e", self.gamma0*(fp*exp(self.kappa_PB*Rep)*(Rep*self.kappa_PB-1)/Rep**2+fm*exp(-self.kappa_PB*Rep)*(-Rep*self.kappa_PB-1)/Rep**2)*ct.elc

    ##mask = self.r*self.lambdaD>res
    ##print mask
    ##plot(self.r[mask]*self.lambdaD, psiNL[mask], 'g',lw=4)
    ### test pour confirmer
    ##self.SetParameters(a=res*1.0e9,I=0.001558,res=400,phi=0.01,sigma=0.056 ) # 0.0366
    ##self.solve()
    ###plot(self.r*self.lambdaD, self.psi, '-g')


    self.ComputeFreeEnergy()
    print("NL Free Energy ", self.A/ct.kT, " kT")
    #self.Plot2D()
    fid = open('cell_a2.dat','w')
    for i in range(len(self.r)):
      fid.write(repr(self.r[i])+' '+repr(self.psi[i])+'\n')
    fid.close()





  def TestPressureExpressions(self):
    close('all')
    figure()
    self = CellModelSimulation()
    i=-1
    Pi_num = []; Pi_out = []; Pi_vir = [];
    phitab = linspace(0.001,0.5,11)
    for phi in phitab:
      i+=1  
      self.SetParameters(a=14.0,I=0.00011526,res=400,phi=phi,sigma=0.007774024696319486,lin=True)
      #self.PrintConstants()
      self.solve()
      self.ComputePressure(); Pi_num.append(self.Pi)
      self.ExactSolutionDHCM(); Pi_out.append(self.Pi)   
      self.ComputeVirialPressure(); Pi_vir.append(self.Pi)
    plot(phitab,Pi_num,'+b',ms=8)
    plot(phitab,Pi_out,'-k',lw=2)
    plot(phitab,Pi_vir,'or',ms=8,mfc='none')
    xlabel('$\\phi$')
    ylabel('$\\Pi_{NI}$')
    i=-1
    Pi_num = []; Pi_out = []; Pi_vir = [];
    phitab = linspace(0.001,0.5,11)
    for phi in phitab:
      i+=1  
      self.SetParameters(a=14.0,I=0.00011526,res=400,phi=phi,sigma=0.003,lin=True)
      #self.PrintConstants()
      self.solve()
      self.ComputePressure(); Pi_num.append(self.Pi)
      self.ExactSolutionDHCM(); Pi_out.append(self.Pi)   
      self.ComputeVirialPressure(); Pi_vir.append(self.Pi)
    plot(phitab,Pi_num,'+b',ms=8)
    plot(phitab,Pi_out,'-k',lw=2)
    plot(phitab,Pi_vir,'or',ms=8,mfc='none')
    i=-1
    Pi_num = []; Pi_out = []; Pi_vir = [];
    phitab = linspace(0.001,0.5,11)
    for phi in phitab:
      i+=1  
      self.SetParameters(a=14.0,I=0.00011526,res=400,phi=phi,sigma=0.02,lin=True)
      #self.PrintConstants()
      self.solve()
      self.ComputePressure(); Pi_num.append(self.Pi)
      self.ExactSolutionDHCM(); Pi_out.append(self.Pi)   
      self.ComputeVirialPressure(); Pi_vir.append(self.Pi)
    plot(phitab,Pi_num,'+b',ms=8)
    plot(phitab,Pi_out,'-k',lw=2)
    plot(phitab,Pi_vir,'or',ms=8,mfc='none')
    xlabel('$\\phi$')
    ylabel('$\\Pi_{NI}$')
    ylim(0,70000)

  
  def TestChargeRegulation(self):
    # Try to reproduce Fig 3 of Gregor Trefalt's review in Langmuir 2015
    # Data for silica     

    close('all')
    self = CellModelSimulation()
    self.SetParameters(a=7.5,I=.01,res=400,phi=.01,pH=8.0,Gamma_tot=array([8.0]),CS=array([2.9]),ksi=array([-1.0]),pK=array([7.5]),lin=False)
    self.PrintConstants()
    self.solve()
    plot(self.r,self.psi)




show()
