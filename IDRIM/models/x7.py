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

PATH = "x7_"

def DeltaTe(gep, Tp, Te, Num, tau_e, Ce):
	'''Electron temperature differential equation.'''
	return (Dt/Ce)*(gep*(Tp-Te)+Num/(2*tau_e))

def DeltaTp(gep, Tp, Te, Num, tau_p, Cp, CON):
	'''Phonon temperature differential equation.'''
	return (Dt/Cp)*((gep*(Te-Tp))+(Num/(2*tau_p))-CON*(1/Tp)*(Tp-T0))
	
def DeltaN(P, Num, tau_e, tau_p):
	'''Non-Thermal energy density differential equation.'''
	return (P-(Num/(2*tau_e))-(Num/(2*tau_p)))*Dt

#Model Action Functions

def Core(I, Relations_dict, SysParams_dict, Fit, CON, time_array=time_array, time_points = time_points):
	'''Model that runs for a given intensity value, relationships, parameters and a Fit constant. Returns dictionary of arrays of values. If you wanted a specific intensity value, could just use this rather than f:Run.'''
	Te=Tp=T0 #Initial Temperature Setting
	Num = 0
	VarArray_dict = GenerateCoreInitial(SysParams_dict)
	
	for k in range(time_points):
		Coef_dict = CoreCoefficientUpdate(Te,Tp, k, I, Relations_dict, SysParams_dict, Fit)
		
		Te = Te + DeltaTe(Coef_dict['gep'], Tp, Te, Num, Coef_dict['tau_e'], Coef_dict['Ce'])
		Tp = Tp + DeltaTp(Coef_dict['gep'], Tp, Te, Num, Coef_dict['tau_p'], Coef_dict['Cp'], CON)
		Num = Num + DeltaN(Coef_dict['S'], Num, Coef_dict['tau_e'], Coef_dict['tau_p'])
		
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
