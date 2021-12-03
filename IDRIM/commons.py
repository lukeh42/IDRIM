#Module Import
from numpy import inf
from scipy.constants import * #constants
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.special import expit
import scipy as sp 
import numpy as np 
from tmm import coh_tmm

from IDRIM.constants import * #IDRIM's constants.py
from IDRIM.relations import *

#Interpolation functions
def Interpolate(Val, X, Y):
	'''Val is value inserted into interpolation function. X, Y are the 2 arrays to interpolate'''
	Interp_Func = interp1d(X, Y)
	return Interp_Func(Val)

def ElectronHC(Temp, Ce):
	'''	Special case of Interpolate function. Requires Electron HC to be generated.	'''
	Interp_Func = interp1d(temperature_array[:-1], Ce) #Special case exists purely due to Ce_array being 1 less dimension than temperature_array
	#could probably merge this and add case detection to f:Interpolate
	return Interp_Func(Temp)

#Transfer Matrix Method
def TMM_Run(n, wavelength, d_list, angle, pol="p"):
	'''Run Transfer Matrix Method. Returns R, T, A in that order.'''
	th0 = angle*degree
	wavelength = wavelength/nm
	n_list = [1, n, 1.5, 1] #List of Refractive Indices
	TMM_out = coh_tmm(pol, n_list, d_list, th0, wavelength)
	R = TMM_out['R']
	T = TMM_out['T']
	A = 1 - R - T
	return R, T, A

def GenInitialRTA(RI_array, wavelength, d_list, angle):
	'''
	Generates initial Reflectance, Transmission and Absorbance coefficients.
	They are also returned in that order. (R, T, A)
	'''
	n0_fresh = Interpolate(300, temperature_array, RI_array)
	nR0, nT0, nA0 = TMM_Run(n0_fresh, wavelength, d_list, angle)
	return nR0, nT0, nA0

def gCoefficient(Tp, Te, wp, pumpf, Fit):
	'''Coupling Coefficient. Returns g_ep. Requires Phonon Temperature.'''
	tau_e = ElectronRelax(300, wp0, pumpf, Fit)#300 comes from Electron Temperature;wp0 is used for plasma freq at 300K
	Lf = tau_e*vf #mean free path calculation
	return (0.562*N*k**2*TD**2*vf)/(Lf*Tp*Ef)

def FitParameterSolver(pumpf, Ce300, Te=300, wp=wp0):
	chi = ((pumpf**2/(4*pi**2*wp))*(1+((2*pi*k*300)/(hbar*pumpf))**2)) #this beats putting this on the next line
	delta = 0.562*((N*k**2*TD**2)/(2*300*Ef*Ce300)) #Tp=300 is where the 300 comes from
	Fit = (chi/gamma0)*(1+delta)
	return Fit

def ElectronRelax(Te, wp, pumpf, Fit):
	'''
	Electron Relaxation Coefficient. Requires Electron Temperature, Plasma Frequency at that temperature.
	'''
	return Fit/((pumpf**2/(4*pi**2*wp))*(1+((2*pi*k*Te)/(hbar*pumpf))**2))

def PhononRelax(g, Ce):
	'''Phonon Relaxation Coefficient. Needs g_ep and Electron HC.'''
	return 2*Ce/g
	
def AbsorbCoeff(n, pumpf):
	'''Absorbance Coefficient. NOT OPTICAL ABSORBANCE! Needs refractive index, or just imag part'''
	return 2*pumpf/c*n.imag
	
def POWER(t, I0, A, alpha, pulse):
	'''Absorbed power density. Needs time point, incident intensity, optical absorption and material absorption'''
	return A*I0*alpha*np.exp(-2*(t/pulse)**2)

def IntensityArrayIndex(I, IntRes):
	'''Returns index for the input intensity value. Intensity arrays must be of correct +1 dimension.'''
	Indexer = int((I/IntRes)) #Works AS LONG AS ARRAYS ARE CREATED WITH +1 DIMENSION!!!
	return Indexer
	
def WavelengthToFrequency(wavelength):
	'''Converts wavelength values to the corresponding frequency value'''
	return (2*pi*c/wavelength)

def GenerateIntensityPoints(IntRes):
	'''Generates number of intensity points given the intensity resolution.'''
	I_points = int(Intensity_max/IntRes)+1
	return I_points

def Normalizer(n, points, Max, Min):
	'''input must be an array'''
	Normal = np.zeros(points)
	for i in range(points):
		Normal[i] = (n[i]-Min)/(Max-Min)
	
	return Normal
