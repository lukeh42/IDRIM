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

'''Model Common Functions'''

def TemporalMatrixAssignment(Dict, Input, j):
	'''Assigns values to dictionary then returns that dictionary given an input dictionary and an index (Intensity index in for loop).'''
	for t in range(time_points):
		Dict['Te'][j][t] = Input['Te'][t]
		Dict['Tp'][j][t] = Input['Tp'][t]
		Dict['Num'][j][t] = Input['Num'][t]
		Dict['T'][j][t] = Input['T'][t]
		Dict['R'][j][t] = Input['R'][t]
		Dict['A'][j][t] = Input['A'][t]
	return Dict
	
def CoreArrayUpdate(Te, Tp, Num, Coefficient_dict, VarArray_dict, i):
	'''Assigns values to dictionary for coefficients and temperatures/Num.'''
	VarArray_dict['Te'][i] = Te
	VarArray_dict['Tp'][i] = Tp
	VarArray_dict['Num'][i] = Num
	VarArray_dict['T'][i] = Coefficient_dict['T']
	VarArray_dict['R'][i] = Coefficient_dict['R']
	VarArray_dict['A'][i] = Coefficient_dict['A']
	return VarArray_dict
	
def CoreCoefficientUpdate(Te, Tp, j, I, Relations_dict, SysParams_dict, Fit):
	'''Function that serves purely to save duplicating this code. Then updates a dictionary with new values, then returns that dictionary.'''
	n = Interpolate(Te, temperature_array, Relations_dict['RI'])
	gep = gCoefficient(Tp)
	Ce_T = ElectronHC(Te, Relations_dict['Ce']) #special case interpolation
	Cp_T = Interpolate(Tp, temperature_array, Relations_dict['Cp']) #old file used DebyePhononHC, slower than this
	wp_T = Interpolate(Te, temperature_array, Relations_dict['wp'])
	tau_e = ElectronRelax(Te, wp_T, WavelengthToFrequency(SysParams_dict['wavelength']), Fit)
	CheckFitParameter(WavelengthToFrequency(SysParams_dict['wavelength']), Fit, ElectronHC(300, Relations_dict['Ce']))#Checking for errors
	tau_p = PhononRelax(gep, Ce_T)
	R, T, A = TMM_Run(n, SysParams_dict['wavelength'], SysParams_dict['thick'], SysParams_dict['angle'])
	CheckART(R, T, A)#checking for errors
	alpha = AbsorbCoeff(n, WavelengthToFrequency(SysParams_dict['wavelength']))
	S = POWER(time_array[j], I*Gwcm2, A, alpha, SysParams_dict['pulse'])
	
	Coef_dict = {'n':n, 'gep':gep,'Ce':Ce_T,'Cp':Cp_T,'wp':wp_T,'tau_e':tau_e,'tau_p':tau_p,'R':R,'T':T,'A':A,'alpha':alpha,'S':S} 
	return Coef_dict
	
def GenerateCoreInitial(Relations_dict, SysParams_dict):
	'''Creates initial arrays and returns them in dictionary. Saves copy paste.'''
	nR0, nT0, nA0 = GenInitialRTA(Relations_dict['RI'], SysParams_dict['wavelength'], SysParams_dict['thick'], SysParams_dict['angle'])
	Te_array = np.full(time_points, T0)
	Tp_array = np.full(time_points, T0)
	Num_array = np.zeros(shape=time_points)
	Tra_array = np.full(time_points, nT0)
	Ref_array = np.full(time_points, nR0)
	Abs_array = np.full(time_points, nA0)
	
	VarArray_dict = {'Te':Te_array, 'Tp':Tp_array, 'Num':Num_array, 'T':Tra_array, 'R':Ref_array, 'A':Abs_array}
	return VarArray_dict

def GenerateMatrixInitial(I_points):
	'''Creates and returns initial matrices. Saves copy pasting it. Then returns them as a dictionary.'''
	Te_matrix = np.zeros(shape=(I_points, time_points))
	Tp_matrix = np.zeros(shape=(I_points, time_points))
	Num_matrix = np.zeros(shape=(I_points, time_points))
	Tra_matrix = np.zeros(shape=(I_points, time_points))
	Ref_matrix = np.zeros(shape=(I_points, time_points))
	Abs_matrix = np.zeros(shape=(I_points, time_points))
	
	VarMatrix_dict = {'Te':Te_matrix, 'Tp':Tp_matrix, 'Num':Num_matrix, 'T':Tra_matrix, 'R':Ref_matrix, 'A':Abs_matrix}
	return VarMatrix_dict
