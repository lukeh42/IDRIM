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
