#Module import
from numpy import inf
from scipy.constants import * #constants
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.special import expit
import scipy as sp 
import numpy as np 

from IDRIM.constants import * #IDRIM's constants.py
from IDRIM.relations import *
from IDRIM.commons import * #Common functions
from IDRIM.checks import * #Checks
from IDRIM.modcom import * #Common model functions

#Differential Equation Definitions

PATH = "x3_"

def DeltaTe(gep, Tp, Te, P, Ce):
	'''Electron temperature differential equation.'''
	return (Dt/Ce)*(gep*(Tp-Te) + P)

def DeltaTp(gep, Tp, Te, Cp):
	'''Phonon temperature differential equation.'''
	return (Dt/Cp)*((gep*(Te-Tp)))
	

#Model Action Functions

def Core(I, Relations_dict, SysParams_dict, Fit, CON, time_array=time_array, time_points = time_points):
	'''Model that runs for a given intensity value, relationships, parameters and a Fit constant. Returns dictionary of arrays of values. If you wanted a specific intensity value, could just use this rather than f:Run.'''
	Te=Tp=T0 #Initial Temperature Setting
	Num = 0 #Rather than rewrite, just eliminate references to Num, so it stays zero.
	VarArray_dict = GenerateCoreInitial(SysParams_dict)
	
	for k in range(time_points):
		Coef_dict = CoreCoefficientUpdate(Te,Tp, k, I, Relations_dict, SysParams_dict, Fit)
		
		Te = Te + DeltaTe(Coef_dict['gep'], Tp, Te, Coef_dict['S'], Coef_dict['Ce'])
		Tp = Tp + DeltaTp(Coef_dict['gep'], Tp, Te, Coef_dict['Cp'])

		
		VarArray_dict = CoreArrayUpdate(Te, Tp, Num, Coef_dict, VarArray_dict, k)
	#end time loop
	return VarArray_dict

def Run(Relations_dict, SysParams_dict, Fit, CON, time_points=time_points):
	'''Model runtime wrapper that cycles through intensities, then outputs desired values to a dictionary of matrices and an array.'''
	#If you wanted to cycle through something other than intensity, probably just modify this wrapper, hell with some adjustments you could probably cycle 1 parameter as well intensity. You'd end up with third order tensors though so good luck have fun.
	I_points = GenerateIntensityPoints(SysParams_dict['Res'])
	Temax_array = np.zeros(I_points)
	VarMatrix_dict = GenerateMatrixInitial(I_points)
	
	for j in range(I_points):
		I = j*SysParams_dict['Res']
		#print(j)
		Hopper = Core(I, Relations_dict, SysParams_dict, Fit, CON) #Hopper is a dictionary as f:Core outputs a dictionary
		
		TemporalMatrixAssignment(VarMatrix_dict, Hopper, j)
		
		Temax_array[j] = max(VarMatrix_dict['Te'][j])
		#print("Temax:", Temax_array)
	#end intensity loop
	VarMatrix_dict.update({'Temax': Temax_array})
	return VarMatrix_dict