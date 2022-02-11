import numpy as np

from IDRIM.constants import * #IDRIM's constants
from IDRIM.relations import *
from IDRIM.commons import *

# Functions around generating relations arrays as well as SMART parameter detection.

def SaveThicknessCase(ParamDict):
	'''Returns dictionary that can be saved. The key part is thick1 and thick2. Can't write to file d_list = [inf, thick1, thick2, inf] so I break it up into the actual parameters.'''
	ParamDict2 = {}
	ParamDict2['wavelength'] = ParamDict['wavelength']
	ParamDict2['pulse'] = ParamDict['pulse']
	ParamDict2['angle'] = ParamDict['angle']
	ParamDict2['Res'] = ParamDict['Res']
	ParamDict2['thick1'] = ParamDict['thick'][1]
	ParamDict2['thick2'] = ParamDict['thick'][2]
	return ParamDict2

def LoadThicknessCase(ParamDict2):
	'''Unused. Basically the reverse of f:SaveThicknessCase. Still here in case it's needed.'''
	ParamDict = {}
	ParamDict['wavelength'] = ParamDict2['wavelength']
	ParamDict['pulse'] = ParamDict2['pulse']
	ParamDict['angle'] = ParamDict2['angle']
	ParamDict['Res'] = ParamDict2['Res']
	ParamDict['thick'][0] = inf
	ParamDict['thick'][1] = ParamDict2['thick1']
	ParamDict['thick'][2] = ParamDict2['thick2']
	ParamDict['thick'][3] = inf
	return ParamDict

def SaveParameters(ParamDictSave):
	'''Saves parameters set in file to history txt file. Uses f:SaveThicknessCase to prevent crash for reasons specified in the aforementioned function.'''
	ParamDict = SaveThicknessCase(ParamDictSave)
	Param_list = (ParamDict['wavelength'], ParamDict['pulse'], ParamDict['angle'], ParamDict['Res'], ParamDict['thick1'], ParamDict['thick2'])
	np.savetxt(Params_path, Param_list) #Params_path is defined in constants. See _documentation.txt for reason
	return

def LoadOldParamsChecker(ParamDictIn):
	'''Function that loads ParameterHistory file, then compares its values with set values and returns 1 if there is no difference, or 0 if there is a difference.'''
	Param_list = np.genfromtxt(Params_path)#Params_path is defined in constants. See _documentation.txt for reason
	ParamDict = SaveThicknessCase(ParamDictIn)
	
	ParamOld = {'wavelength':Param_list[0], 'pulse':Param_list[1], 'angle':Param_list[2], 'Res':Param_list[3], 'thick1':Param_list[4], 'thick2':Param_list[5]}

	if ParamOld == ParamDict: #directly compares the dictionaries. They must be identical.
		Outline = 1
	else:
		Outline = 0
		
	return Outline
	
def Regeneration(SysParams_dict, OVERRIDE = 0):
	'''Master function that calls a bunch of other functions dependent on outcome of f:LoadOldParamsChecker. Optional OVERRIDE is avaliable to force a regeneration of arrays regardless. Set to 1 (or anything that isnt 0)
	Override = -1 skips loading, saving and RI generation. Used only for testing. '''
	if OVERRIDE != -1:
		Skipper = LoadOldParamsChecker(SysParams_dict)
	if OVERRIDE != 0:
		Skipper = 2
	if Skipper == 1:
		print("No change in parameters detected.")
	elif Skipper == 0:
		print("Change in parameters detected. Regenerating arrays.")
	else:
		print("OVERRIDE Detected. Regenerating.")
		Skipper = 0
	mu_array = GenMu(Skipper) #strictly speaking Mu only relies on max temp and DIM, optimisation?
	Cp_array = GenPhononHC(Skipper)
	Ce_array = GenElectronHC(mu_array, Skipper)
	wp_array = GenPlasmaFrequency(GenAvgEffMass(mu_array), Skipper)
	if OVERRIDE != -1:
		RI_array = GenRefractiveIndex(WavelengthToFrequency(SysParams_dict['wavelength']), wp_array, Skip=Skipper)
	if OVERRIDE != -1:
		SaveParameters(SysParams_dict)
		Ce300 = ElectronHC(300, Ce_array)
		Fit = FitParameterSolver(WavelengthToFrequency(SysParams_dict['wavelength']), Ce300)
	if OVERRIDE == -1:
		RI_array = 0
		Fit = 0
	
	print("\n\nREGENERATION COMPLETE.")
	
	return mu_array, Cp_array, Ce_array, wp_array, RI_array, Fit
	
