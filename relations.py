#Module imports
from numpy import inf
from scipy.constants import * #constants
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.special import expit
import scipy as sp 
import numpy as np 
from IDRIM.constants import * #IDRIM's constants.py

#Dispersion Relations

def TCODispersion(E):
	'''Returns k value of TCO Dispersion relation for an input energy value.'''
	return np.sqrt((2*me_min/hbar**2)*(E+eta*E**2))

def ParabolicDispersion(E):
	'''Returns k value of parabolic dispersion relation for an input energy value.'''
	return np.sqrt(2*me_min*E/hbar**2)

def GenDispersionArrays():
	'''Generates and returns dispersion arrays using GLOBAL energy array.'''
	TCO_k_array = TCODispersion(energy_array)
	Parabolic_k_array = ParabolicDispersion(energy_array)
	return TCO_k_array, Parabolic_k_array

#Density of States

def TCO_DOS(E):
	'''Returns density of states for input energy value for TCO.'''
	return (1/(2*pi**2))*(2*me_min/hbar**2)**1.5*(E+eta*E**2)**0.5*(1+2*eta*E)

def ParabolicDOS(E):
	'''Returns density of states for an input energy for parabolic case.'''
	return (1/(2*pi**2))*(2*me_min/hbar**2)**1.5*E**0.5

def GenDOSArrays():
	'''Generates and returns TCO/Parabolic DOS arrays in that order using GLOBAL energy array.'''
	TCO_DOS_array = TCO_DOS(energy_array)
	Parabolic_DOS_array = ParabolicDOS(energy_array)
	return TCO_DOS_array, Parabolic_DOS_array
	
#Fermi Dirac Distribution

def FermiDirac(E, mu_T, T):
	'''Fermi Dirac Distribution. Returns value given input parameters.'''
	x = (mu_T - E)/(k*T)
	return expit(x) #expit(x)# 1/(1+exp(-x))

#Chemical Potential

def Mu_Integrand(E, mu_T, T):
	'''Returns integrand of integral for calculating chemical potential. You probably want f:Mu_Integrate.'''
	return FermiDirac(E, mu_T, T)*TCO_DOS(E)

def Mu_Integrate(mu_T, T):
	'''Returns integration result of chemical potential.'''
	return quad(Mu_Integrand, 0, 100*e, args=(mu_T, T))[0]

def Mu_Function(mu_T, T):
	'''Returns result of chemical potential integral minus the electron density constant. You probably want f:GenMu.'''
	return Mu_Integrate(mu_T, T)-N
	

def GenMu(Skip=0):
	'''Generates and returns chemical potential. Optional Skip arg controls whether to Load from saved file (Skip=1) or Generate Mu.(Skip=else)'''
	Mu_Path = "data\Mu_array.txt"
	if Skip == 1:
		print("Not solving for Chemical Potential. Loading from previously saved data.")
		mu_array = np.genfromtxt(Mu_Path)
		print("Mu Loaded.")
	
	else:
		print("Solving Chemical Potential. This may take some time.")
		mu_array = np.zeros(shape=DIM)
		for i in range(DIM):
			T = temperature_array[i]
			mu_array[i] = fsolve(Mu_Function, Ef, args=(T))
			#mu_array = np.append(mu_array, fsolve(Mu_Function, 0.5*e, args=(T)))
			
		np.savetxt(Mu_Path, mu_array)
		print("Mu Generated and Saved.")
	return mu_array
	
# Phonon Heat Capacity

def DebyeIntegrand(x, T):
	'''Defines the integrand in the Debye model. You probably want f:DebyePhononHC.'''
	return 3*9.0*N*k*((T/TD)**3)*((x**4 *np.exp(x))/((np.exp(x)-1)**2))
	#don't know why the 3 is in there.
def DebyePhononHC(T):
	'''Returns integral result of Debye Model for Phonon Heat Capacity.'''
	return quad(DebyeIntegrand, 0, (TD/T), args=(T))[0]

def GenPhononHC(Skip=0):
	'''Generates Phonon Heat Capacity using the Debye Model. Skip arg controls whether it is loaded from file (Skip=1) or generated from scratch (Skip=else).'''
	Cp_path = "data\Cp_array.txt"
	if Skip == 1:
		print("Not solving Phonon Heat Capacity. Loading from saved file.")
		Cp_array = np.genfromtxt(Cp_path)
		print("Phonon Heat Capacity loaded.")
	else:
		print("Generating Phonon Heat Capacity.")
		Cp_array = np.zeros(shape=DIM)
		for i in range(DIM):
			T = temperature_array[i]
			Cp_array[i] = DebyePhononHC(T)
		np.savetxt(Cp_path, Cp_array)
		print("Phonon Heat Capacity generated and saved.")
	return Cp_array
	
# Electron Heat Capacity

def ElectronHCIntegrand(E, mu_T, T):
	'''Integrand for calculating energy. You probably want f:ElectronHCEnergy.'''
	return TCO_DOS(E)*E*FermiDirac(E, mu_T, T)
	
def ElectronHCEnergy(E, mu_T, T):
	'''Returns integration result for use in calculating electron heat capacity.'''
	return quad(ElectronHCIntegrand, 0, 100*e, args=(mu_T, T))[0]

def GenElectronHC(mu_array, Skip=0):
	'''Generates electron heat capacity. Requires chemical potential input. Skip=1 loads a pre-generated file and Skip=0 generates from scratch.'''
	Ce_path = "data\Ce_array.txt"
	if Skip == 1:
		print("Not solving Electron Heat Capacity. Loading from saved file.")
		Ce_array = np.genfromtxt(Ce_path)
		print("Electron Heat Capacity loaded.")
	else:
		print("Generating Electron Heat Capacity.")
		HC_energy_array = np.zeros(shape=DIM)
		for i in range(DIM):
			mu_T = mu_array[i]
			T = temperature_array[i]
			HC_energy_array[i] = ElectronHCEnergy(energy_array, mu_T, T)
		Ce_array = ((np.diff(HC_energy_array))/(np.diff(temperature_array)))
		np.savetxt(Ce_path, Ce_array)
		print("Electron Heat Capacity generated and saved.")
	return Ce_array

# Effective Masses

def EffMassExtended(E):
	'''Returns Standard effective mass'''
	return me_min*(1+2*eta*E)

def EffMassOptical(E):
	'''Returns Optical effective mass'''
	return (a**2/(hbar**2*b*(E+a)))**(-1)

def EffMassConventional(E):
	'''Returns conventional effective mass'''
	return ((a**4)/(hbar**2*b*((E+a))**3) )**(-1)

def EffMassWeighted(E):
	'''Weighted effective mass. 2/3 optical + 1/3 conventonal. Used in this model.'''
	return (((2/3)*a**2/(hbar**2*b*(E+a))) + ((1/3)*(a**4)/(hbar**2*b*((E+a))**3)))**(-1)

def GenEffMasses():
	'''Generates effective masses and returns their arrays. ORDER: EXT, OPT, CON, WEI.'''
	Extended_array = EffMassExtended(energy_array)
	Optical_array = EffMassOptical(energy_array)
	Conventional_array = EffMassConventional(energy_array)
	Weighted_array = EffMassWeighted(energy_array)
	return Extended_array, Optical_array, Conventional_array, Weighted_array

#Average Effective Mass

def AvgEffMassIntegrand(E, mu_T, T):
	'''Average Effective Mass integrand. You probably want f:AvgEffMassIntegrate.'''
	return (1/N)*(TCO_DOS(E)*FermiDirac(E, mu_T, T))/(EffMassWeighted(E))

def AvgEffMassIntegrate(mu_T, T):
	'''Integrates for average effective mass.'''
	return 1/(quad(AvgEffMassIntegrand, 0, 100*e, args=(mu_T, T))[0])

def GenAvgEffMass(mu_array):
	'''Generates and returns average effective mass. Requires chemical potential input.'''
	AvgEffMass_array = np.zeros(shape=DIM)
	for i in range(DIM):
		mu_T = mu_array[i]
		T = temperature_array[i]
		AvgEffMass_array[i] = AvgEffMassIntegrate(mu_T, T)
	return AvgEffMass_array
	
# Plasma Frequency

def GenPlasmaFrequency(AEM, Skip=0):
	'''Generates and returns plasma frequency. Requires average effective mass array. Skip=0 generates, Skip=1 loads from saved file.'''
	wp_path = "data\wp_array.txt"
	if Skip == 1:
		print("Loading Plasma Frequency from saved file.")
		wp_array = np.genfromtxt(wp_path)
		print("Loaded Plasma Frequency.")
	else:
		print("Generating Plasma Frequency")
		wp_array = np.sqrt((N*e**2)/(epsilon_0*AEM))	
		np.savetxt(wp_path, wp_array)
		print("Plasma Frequency array generated and saved.")
	return wp_array

# Refractive Index

def Permittivity(w, wp, gamma=gamma0):
	'''Simply returns permittivity given a frequency, plasma frequency and optionally a scattering rate.'''
	return eps_inf - ((wp**2)/(w**2+1j*w*gamma))

def GenRefractiveIndex(pump_freq, wp, gamma=gamma0, Skip=0):
	'''Generates refractive index array. Requires pump frequency, plasma frequency array. Skip=0 generates from scratch, Skip=1 loads from data.'''
	#Refractive Index must be split into real and imaginary when saving to file. np.savetxt saves complex numbers as one object, ie "(REAL + IMAGj)" where the brackets ARE included. np.genfromtxt then creates a bunch of strings or whatever.
	#Hence the only way to save and load refractive index is to split up real and imaginary parts and then recombine them after loading both files. Kinda interesting bug actually
	RIr_path = "data\RIr_array.txt"
	RIi_path = "data\RIi_array.txt"
	if Skip == 1:
		print("Loading Refractive Index from saved file.")
		RIr_array = np.genfromtxt(RIr_path)
		RIi_array = np.genfromtxt(RIi_path)
		RI_array = RIr_array + 1j*RIi_array
		print("Loaded Refractive Index")
	else:
		print("Generating Refractive Index.")
		RI_array = np.sqrt(Permittivity(pump_freq, wp, gamma))
		np.savetxt(RIr_path, RI_array.real)
		np.savetxt(RIi_path, RI_array.imag)
		print("Refractive Index array generated and saved.")
	return RI_array


