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
from IDRIM.ymodcom import * #ymodel common funcs

PATH = "y1_"
CON = 1e18

def DeltaN(N, P, tep, tee, ycone, yconp):
	return (P - ((N)/(yconp*tep)) - ((N)/(ycone*tee)))*Dt

def DeltaTe(Ce, gep, tee, N, Te, Tp, ycone):
	return (gep*(Tp-Te) + ((N)/(ycone*tee)))*(Dt/Ce)

def DeltaTp(Cp, gep, tep, N, Te, Tp, yconp):
	return (gep*(Te-Tp) + ((N)/(yconp*tep)))*(Dt/Cp)

def DeltaTpAlt(Cp, gep, tep, N, Te, Tp, yconp):
	return (gep*(Te-Tp) + ((N)/(yconp*tep))-CON*(Tp-300))*(Dt/Cp)


def Core(I, Relations, Parameters, Y):
	Te=Tp=T0
	Num = 0
	VarArray_dict = GenerateCoreInitial(Relations, Parameters)
	for t in range(time_points):
		
		wp_T = Interpolate(Te, temperature_array, Relations['wp'])#needed for Y13 so done here
		
		if Y['Y13'] == 1: #Model Y-13: Non-Constant Drude Scattering Rate
			gamma = ScatteringTemp(Te) #has unusued second arg, defaulting to 15,000K (as per EM paper)
			n = np.sqrt(Permittivity(WavelengthToFrequency(Parameters['wavelength']), wp_T, gamma))
		else:
			n = Interpolate(Te, temperature_array, Relations['RI']) #Relations['RI'] will be constant gamma
		
		gep = gCoefficient(Tp, WavelengthToFrequency(Parameters['wavelength']), Relations['Fit'])
		
		R, T, A = TMM_Run(n, Parameters['wavelength'], Parameters['thick'], Parameters['angle'])
		alpha = AbsorbCoeff(n, WavelengthToFrequency(Parameters['wavelength']))
		tee = ElectronRelax(Te, wp_T, WavelengthToFrequency(Parameters['wavelength']), Relations['Fit'])
		
		if Y['Y2'] == 1: #Model Y-2: Simple Electron Heat Capacity
			Ce = yElectronHC(Te)
		else:
			Ce = ElectronHC(Te, Relations['Ce'])
		
		if Y['Y3'] == 1: #Model Y-3: Constant Phonon Heat Capacity
			Cp = yCp
		else:
			Cp = Interpolate(Tp, temperature_array, Relations['Cp'])
		
		if Y['Y5'] == 1: #Model Y-5: Power Function Change
			Ip = 20.8*Gwcm2
			alphaprime=4.18e6
			P = yPower(time_array[t], I*Gwcm2, alphaprime, Parameters['pulse'])
		else:
			P = POWER(time_array[t], I*Gwcm2, A, alpha, Parameters['pulse'])
		
		if Y['Y7'] == 1: #Model Y-7: Phonon Relaxation Coefficient Change
			tep = yPhononRelax(Ce, gep)
			yconp = 1
			ycone = 1 #yconp and ycone are separate as might want to have them be separate in the future
		else:
			tep = PhononRelax(gep, Ce)
			yconp = 2
			ycone = 2
		
		#Model Y-11 occurs further down, it's phonon decay term.
		
		#CheckFitParameter(WavelengthToFrequency(Parameters['wavelength']), Relations['Fit'], ElectronHC(300, Relations['Ce']))#Checking for errors
		CheckART(R, T, A) # checking for errors
		#print(Relations['Fit'])
		Te = Te + DeltaTe(Ce, gep, tee, Num, Te, Tp, ycone)
		if Y['Y11'] == 1:
			Tp = Tp + DeltaTpAlt(Cp, gep, tep, Num, Te, Tp, yconp)
		else:
			Tp = Tp + DeltaTp(Cp, gep, tep, Num, Te, Tp, yconp)
		Num = Num + DeltaN(Num, P, tep, tee, ycone, yconp)

		Coef_dict = {'n':n, 'gep':gep,'Ce':Ce,'Cp':Cp,'wp':wp_T,'tau_e':tee,'tau_p':tep,'R':R,'T':T,'A':A,'alpha':alpha,'S':P} 
		VarArray_dict = CoreArrayUpdate(Te, Tp, Num, Coef_dict, VarArray_dict, t)
	return VarArray_dict
	
	
def Intensity_Run(Relations, Parameters, Y):
	I_points = GenerateIntensityPoints(Parameters['Res'])
	Temax_array = np.zeros(I_points)
	VarMatrix_dict = GenerateMatrixInitial(I_points)
	
	for j in range(I_points):
		I = j*Parameters['Res']
		Hopper = Core(I, Relations, Parameters, Y)
		
		TemporalMatrixAssignment(VarMatrix_dict, Hopper, j)
		
		Temax_array[j] = max(VarMatrix_dict['Te'][j])
	#ending intensity loop
	VarMatrix_dict.update({'Temax': Temax_array})
	return VarMatrix_dict
	
def Angle_Run(Relations, Parameters, I, Y, AngleRes=10):
	Angle_max = 90
	Angle_points = int(Angle_max/AngleRes)
	VarMatrix_dict = GenerateMatrixInitial(Angle_points)
	Temax_array = np.zeros(Angle_points)
	
	for o in range(Angle_points):
		O = o*AngleRes
		Parameters.update({'angle': O})
		print("Current Angle:", O)
		Hopper = Core(I, Relations, Parameters, Y)
		
		TemporalMatrixAssignment(VarMatrix_dict, Hopper, o)
		
		Temax_array[o] = max(VarMatrix_dict['Te'][o])
		
		#end angle loop
	VarMatrix_dict.update({'Temax': Temax_array})
	return VarMatrix_dict
	
def Wavelength_Run(Relations, Parameters, I, Y, WavelengthRes=100):
	Wavelength_min = 700
	Wavelength_max = 1500
	#enz wavelength depends on sample
	
	Lambda_points = int((Wavelength_max - Wavelength_min)/WavelengthRes)
	VarMatrix_dict = GenerateMatrixInitial(Lambda_points)
	Temax_array = np.zeros(Lambda_points)
	
	for l in range(Lambda_points):
		L = (Wavelength_min + l*WavelengthRes)*nm
		Parameters.update({'wavelength': L})
		Relations.update({'RI': GenRefractiveIndex(WavelengthToFrequency(Parameters['wavelength']), Relations['wp'], Skip=0)})
		Hopper = Core(I, Relations, Parameters, Y)
		
		TemporalMatrixAssignment(VarMatrix_dict, Hopper, l)
		
		Temax_array[l] = max(VarMatrix_dict['Te'][l])
		#end loop
	VarMatrix_dict.update({'Temax': Temax_array})
	return VarMatrix_dict

def Pulse_Run(Relations, Parameters, I, Y, PulseRes=5):
	Pulse_min = 90
	Pulse_max = 120
	Pulse_points = int((Pulse_max - Pulse_min)/PulseRes)
	VarMatrix_dict = GenerateMatrixInitial(Pulse_points)
	
	for p in range(Pulse_points):
		P = p*PulseRes*1e-15
		Parameters.update({'pulse': P})
		Hopper = Core(I, Relations, Parameters, Y)
		
		TemporalMatrixAssignment(VarArray_dict, Hopper, p)
		
		Temax_array[p] = max(VarMatrix_dict['Te'][p])
		#end loop
	VarMatrix_dict.update({'Temax': Temax_array})
	return VarMatrix_dict
