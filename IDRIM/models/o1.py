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

#Different definition of Electron Heat Capacity
def oElectronHC(Te):
	gamma = (3*pi**2*7.31e26*k)/(np.sqrt(36*Tf**2+4*pi**4*Te**2))
	Ce = gamma*Te
	return Ce

def oPhononRelax(g, Ce):
	return Ce/g

def oPower(t, pulse, Ip=20.8, alpha=4.18e6):
	return Ip*Gwcm2*alpha*np.exp(-2*(t/pulse)**2)
#Differential Equation Definitions

PATH = "o1_"

def DeltaTe(gep, Tp, Te, Num, tau_e, Ce):
	'''Electron temperature differential equation.'''
	return (Dt/Ce)*(gep*(Tp-Te)+Num/(2*tau_e))

def DeltaTp(gep, Tp, Te, Num, tau_p, Cp):
	'''Phonon temperature differential equation.'''
	return (Dt/Cp)*((gep*(Te-Tp))+(Num/(2*tau_p)))
	
def DeltaN(P, Num, tau_e, tau_p):
	'''Non-Thermal energy density differential equation.'''
	return (P-(Num/(2*tau_e))-(Num/(2*tau_p)))*Dt

#Model Action Functions

def Core(I, Relations_dict, SysParams_dict, Fit, time_array=time_array, time_points = time_points):
	'''Model that runs for a given intensity value, relationships, parameters and a Fit constant. Returns dictionary of arrays of values. If you wanted a specific intensity value, could just use this rather than f:Run.'''
	Te=Tp=T0 #Initial Temperature Setting
	Num = 0
	VarArray_dict = GenerateCoreInitial(SysParams_dict)
	
	for k in range(time_points):
	
		n = Interpolate(Te, temperature_array, Relations_dict['RI'])
		gep = gCoefficient(Tp)
		Ce_T = oElectronHC(Te) #new electron heat capacity
		Cp_T = 2.6e6 #Constant
		wp_T = Interpolate(Te, temperature_array, Relations_dict['wp'])
		tau_e = ElectronRelax(Te, wp_T, WavelengthToFrequency(SysParams_dict['wavelength']), Fit)
		CheckFitParameter(WavelengthToFrequency(SysParams_dict['wavelength']), Fit, oElectronHC(300))#Checking for errors
		tau_p = oPhononRelax(gep, Ce_T)
		R, T, A = TMM_Run(n, SysParams_dict['wavelength'], SysParams_dict['thick'], SysParams_dict['angle'])
		CheckART(R, T, A)#checking for errors
		alpha = AbsorbCoeff(n, WavelengthToFrequency(SysParams_dict['wavelength']))
		S = oPower(k, pulse=SysParams_dict['pulse'])
	
		Coef_dict = {'n':n, 'gep':gep,'Ce':Ce_T,'Cp':Cp_T,'wp':wp_T,'tau_e':tau_e,'tau_p':tau_p,'R':R,'T':T,'A':A,'alpha':alpha,'S':S} 
		
		Te = Te + DeltaTe(Coef_dict['gep'], Tp, Te, Num, Coef_dict['tau_e'], Coef_dict['Ce'])
		Tp = Tp + DeltaTp(Coef_dict['gep'], Tp, Te, Num, Coef_dict['tau_p'], Coef_dict['Cp'])
		Num = Num + DeltaN(Coef_dict['S'], Num, Coef_dict['tau_e'], Coef_dict['tau_p'])
		
		VarArray_dict = CoreArrayUpdate(Te, Tp, Num, Coef_dict, VarArray_dict, k)
	#end time loop
	return VarArray_dict

def Run(Relations_dict, SysParams_dict, time_points=time_points):
	'''Model runtime wrapper that cycles through intensities, then outputs desired values to a dictionary of matrices and an array.'''
	#If you wanted to cycle through something other than intensity, probably just modify this wrapper, hell with some adjustments you could probably cycle 1 parameter as well intensity. You'd end up with third order tensors though so good luck have fun.
	I_points = GenerateIntensityPoints(SysParams_dict['Res'])
	Temax_array = np.zeros(I_points)
	VarMatrix_dict = GenerateMatrixInitial(I_points)
	
	Ce_new = oElectronHC(temperature_array)
	Fit_new = FitParameterSolver(WavelengthToFrequency(SysParams_dict['wavelength']), oElectronHC(300))
	
	Relations_dict.update({'Ce': Ce_new})
	Relations_dict.update({'Fit': Fit_new})
	
	for j in range(I_points):
		I = j*SysParams_dict['Res']
		#print(j)
		Hopper = Core(I, Relations_dict, SysParams_dict, Fit_new) #Hopper is a dictionary as f:Core outputs a dictionary
		
		TemporalMatrixAssignment(VarMatrix_dict, Hopper, j)
		
		Temax_array[j] = max(VarMatrix_dict['Te'][j])
		#print("Temax:", Temax_array)
	#end intensity loop
	VarMatrix_dict.update({'Temax': Temax_array})
	return VarMatrix_dict
