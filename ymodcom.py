from numpy import inf
from scipy.constants import * #constants
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.special import expit
import scipy as sp 
import numpy as np 

from IDRIM.constants import *
from IDRIM.relations import *
from IDRIM.checks import *
from IDRIM.commons import *

def yElectronHC(Te): # Model Y-2
	return Te*((3*pi**2*N*k)/(np.sqrt(36*Tf**2 + 4*pi**4*Te**2)))

yCp = 2.6e6 # Model Y-3

def yPower(t, Ip, alpha, pulse): # Model Y-5
	return Ip*alpha*np.exp(-2*(t/pulse)**2)

def yPhononRelax(Ce, gep): # Model Y-7
	return Ce/gep

def ScatteringTemp(Te, gaT0=15000):
'''Temperature dependent drude scattering/damping coefficient. Default for T0 from EM paper [Time Domain ENZ]'''
	return gamma0*(1+Te/gaT0)

'''
Depreciated functions for calculating plasma frequency through another method. 
After finding that it produces the same result as the current method, nothing more was done.
They exist here as legacy code, in case they are needed again.
'''
def FermiDiracDifferential(E, mu, T):
	return (1/(k*T))*(np.exp((E-mu)/(k*T)))/(((np.exp((E-mu)/(k*T))+1))**2)

def EMPlasmaFrequencyIntegrand(E, mu, T):
	return ((e**2)/(3*me_min*epsilon_0*pi**2))*((2*me_min)/(hbar**2))**(1.5)*((E+eta*E**2)**(1.5))*((1+2*eta*E)**(-1))*FermiDiracDifferential(E, mu, T)
	
def EMPlasmaFrequencyIntegrate(mu, T):
	return quad(EMPlasmaFrequencyIntegrand, 0, 100*e, args=(mu, T))[0]

def EMPlasmaFrequencyGenerator(mu_array):
	'''Generates and returns average effective mass. Requires chemical potential input.'''
	EM_Plasma_array = np.zeros(shape=DIM)
	for i in range(DIM):
		mu = mu_array[i]
		T = temperature_array[i]
		EM_Plasma_array[i] = np.sqrt(EMPlasmaFrequencyIntegrate(mu, T))
	return EM_Plasma_array
	
