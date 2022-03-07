import numpy as np
import scipy as sp
from scipy.special import expit
from scipy.integrate import quad
from scipy.constants import *

from IDRIM.constants import *
from IDRIM.relations import *

import matplotlib.pyplot as plt
'''Different Mu functions, contains N as variable.'''
def loc_Mu_Integrand(E, mu_T, T):
	return FermiDirac(E, mu_T, T)*TCO_DOS(E)

def loc_Mu_Integrate(mu_T, T):
	return quad(Mu_Integrand, 0, 100*e, args=(mu_T, T))[0]

def loc_Mu_Function(mu_T, T, N):
	return Mu_Integrate(mu_T, T)-N


'''N and m0 funcs'''
def NumberDensityIntegrand(E, m0):
	return (1/(2*pi**2))*(2*m0/hbar**2)**1.5*(E+eta*E**2)**0.5*(1+2*eta*E)

def CalculateNumberDensity(m0, Ef=1*e): #!M!
	'''Input minimum effective mass to calculate number density at 0 kelvin. Uses constant Fermi Energy (1eV)'''
	return quad(NumberDensityIntegrand, 0, Ef, args=m0)[0]
	
def ModdedWeightedEffMass(E):
	return (((2*a**2)/(3*hbar**2*(E+a)))+(a**4)/(3*hbar**2*(E+a)**3))
	
def IntegrandC1a(E, mu, T):
	return ((E+eta*E**2)**0.5)*(1+2*eta*E)*FermiDirac(E, mu, T)*ModdedWeightedEffMass(E)

def IntegrateC1a(mu, T):
	return quad(IntegrandC1a, 0, 100*e, args=(mu, T))[0]

def IntegrandC2(E, mu, T):
	return ((E+eta*E**2)**0.5)*(1+2*eta*E)*FermiDirac(E, mu, T)

def IntegrateC2(mu, T):
	return quad(IntegrandC2, 0, 100*e, args=(mu, T))[0]

def CalculateMass(wp0, Ef=1*e): #!M!
	'''Calculates minimum effective mass given plasma frequency at 300K'''
	C2 = IntegrateC2(Ef, 300)
	C1 = (IntegrateC1a(Ef, 300))/(C2)
	return wp0**4 * ((epsilon_0**2 *pi**4)/(C1**2 * C2**2 * e**4 * eta**2 *hbar**4))*((hbar**2)/(2))**3
	
'''Effective Mass Functions'''
def loc_EffMassWeighted(E, b):
	return (((2/3)*a**2/(hbar**2*b*(E+a))) + ((1/3)*(a**4)/(hbar**2*b*((E+a))**3)))**(-1)

def loc_AvgEffMassIntegrand(E, mu_T, T, b):
	return (1/N)*(TCO_DOS(E)*FermiDirac(E, mu_T, T))/(loc_EffMassWeighted(E, b))

def loc_AvgEffMassIntegrate(mu_T, T, b):
	return 1/(quad(loc_AvgEffMassIntegrand, 0, 100*e, args=(mu_T, T, b))[0])

def ENZFreqSolver(wp, epsilon_infinity, gamma):
	return ((1/(2*pi))*np.sqrt((wp**2)/(epsilon_infinity) - gamma**2))

'''Solution and Graphing'''
def Solver(wp0, Ef=1*e, ENZ=212e12):
	'''Solves everything given plasma frequency at 300K. Supports custom Fermi energy and ENZ frequency'''
	Frequency_array = np.linspace(150, 300, 1000)*1e12*2*pi
	PermExperiment = Permittivity(Frequency_array, wp0)
	
	m0 = CalculateMass(wp0, Ef)
	
	bloc =  (m0)/(2*eta*hbar**2)
	
	print("Minimum Effective Mass:", m0/m_e, "m_e")
	Nlocal = CalculateNumberDensity(m0, Ef)
	print("Number Density at 0 Kelvin:",Nlocal)
	mu300 = fsolve(loc_Mu_Function, 0.5*e, args=(300, Nlocal))
	print("Chemical Potential at 300K:", mu300/eV, "eV")
	
	m_avg = loc_AvgEffMassIntegrate(mu300, 300, bloc) #Calculate average effective mass.
	PlasmaFreq = np.sqrt((N*e**2)/(epsilon_0*m_avg))
	PermModel = Permittivity(Frequency_array, PlasmaFreq)
   
	plt.plot(Frequency_array/(1e12*2*pi), PermExperiment.real, label="300K Plasma Freq Constant")
	plt.plot(Frequency_array/(1e12*2*pi), PermModel.real, label="Calculated Plasma Frequency")
	plt.axhline(0, color="black", linestyle="dashed", label="ENZ Line")
	plt.axvline(ENZ/1e12, color="black", linestyle="dotted", label="ENZ Frequency")
	plt.xlabel("Frequency (THz)")
	plt.ylabel("Real(Permittivity)")
	plt.title('Model Fermi Energy: %f eV' % (Ef/eV))
	plt.legend()
	return