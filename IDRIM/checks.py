#Module Imports
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
from IDRIM.commons import *

def CheckART(R, T, A):
	'''Checks that optical coefficients equal unity. This function is ran a lot inside the model. Will raise a ValueError if something goes wrong.'''
	Unity = round(R+T+A, 3) #rounds optical coefficient sum to 3sf so it'll be 1 if it's ever so slightly over due to rounding when RTA are generated
	if Unity != 1:
		raise ValueError('Optical Coefficients are not summing to 1. Something very bad has occured.', R, T, A)
	return

def CheckFitParameter(pumpf, Fit, Ce300):
	'''Checks that the Fit constant is set such that it's defining equation is unity. Raises a ValueError if this is not the case.'''
	dum_Unity = ((1/ElectronRelax(300, wp0, pumpf, Fit)) + (1/PhononRelax(gCoefficient(300), Ce300)))/gamma0
	Unity = round(dum_Unity, 1) #rounds for same reason as CheckART.
	if Unity != 1:
		raise ValueError('Fit parameter is not being set correctly.', 'Unity:',dum_Unity, 'Fit:',Fit)
	return

def CheckRefractiveIndex(RI_array):
	'''Checks refractive index, should be approximately the constant RI. Requires RI_array'''
	check_n0 = Interpolate(300, temperature_array, RI_array)
	print("Constant n0:", n0)
	print("Generated n0:", check_n0)
	return 

def OldMatrixMethod(n_ito, n_sub, freq, d_ito, d_sub, points):
	'''
	n_ito: RI array of ITO sample. n_sub: RI of substrate. freq: frequency array. d_ito: thickness of ito sample.
	d_sub: thickness of substrate. points: length of arrays.
	This method does not account for angle, so when comparing to TMM module, set the angle to 0,
	or use f:MatrixMethod.
	'''
	def TrigInsider(RI, f, thick):#only used here so local func definition
		return RI*f*thick/c
	
	T = np.zeros(points)
	R = np.zeros(points)
	A = np.zeros(points)
	
	for i in range(points):
		in_ito = TrigInsider(n_ito[i], freq[i], d_ito)
		Z_ito = np.array([[np.cos(in_ito), -(1j/n_ito[i])*np.sin(in_ito)],[-1j*n_ito[i]*np.sin(in_ito), np.cos(in_ito)]])
		
		in_sub = TrigInsider(n_sub, freq[i], d_sub)
		Z_sub = np.array([[np.cos(in_sub),(-1j/n_sub)*np.sin(in_sub)],[-1j*n_sub*np.sin(in_sub), np.cos(in_sub)]])
		
		Z = np.matmul(Z_sub, Z_ito)
		
		r = (Z[0][0] + Z[1][0] - Z[0][1] - Z[1][1])/(Z[0][0] + Z[0][1] + Z[1][0] + Z[1][1])
		t = (2)/(Z[0][0] + Z[0][1] + Z[1][0] + Z[1][1])
		
		T[i] = np.abs(t)**2
		R[i] = np.abs(r)**2
		A[i] = 1 - R[i] - T[i]
		
	return R, T, A

def MatrixMethod(n_ito, n_sub, wavel, d_ito, d_sub, points):
	'''
	n_ito: RI array of ITO sample. n_sub: RI of substrate. freq: frequency array. d_ito: thickness of ito sample.
	d_sub: thickness of substrate. points: length of arrays.
	Use this in conjunction with f:OldMatrixMethod to compare the two.
	'''
	T = np.zeros(points)
	R = np.zeros(points)
	A = np.zeros(points)	
	
	for i in range(points):
		R[i], T[i], A[i] = TMM_Run(n_ito[i], wavel[i], [inf, d_ito/nm, d_sub/nm, inf], 0)
	return R, T, A