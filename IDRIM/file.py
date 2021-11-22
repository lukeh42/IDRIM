import numpy as np
from IDRIM.constants import *

#Model output saving and loading.
def ModelOutputSave(PATH, M_dict):
	'''Saves a model output dictionary to files specified in PATH parameter. Note that the PATH parameter can also be an expression. So eg, ModelOutputSave(x1.PATH+"trial41591", dictionary) is totally valid and will work! No override checks. Be careful.'''
	print("Saving files in", PATH, "directory.")
	np.savetxt(PATH+'_Te.txt', M_dict['Te'])
	np.savetxt(PATH+'_Tp.txt', M_dict['Tp'])
	np.savetxt(PATH+'_Num.txt', M_dict['Num'])
	np.savetxt(PATH+'_T.txt', M_dict['T'])
	np.savetxt(PATH+'_R.txt', M_dict['R'])
	np.savetxt(PATH+'_A.txt', M_dict['A'])
	return
	
def ModelOutputLoad(PATH):
	'''Loads a model output dictionary from files specified in PATH parameter. No checks if it is actually there. Be careful.'''
	print("Loading files from", PATH, "directory.")
	M_dict = {}
	M_dict['Te'] = np.genfromtxt(PATH+'_Te.txt')
	M_dict['Tp'] = np.genfromtxt(PATH+'_Tp.txt')
	M_dict['Num'] = np.genfromtxt(PATH+'_Num.txt')
	M_dict['T'] = np.genfromtxt(PATH+'_T.txt')
	M_dict['R'] = np.genfromtxt(PATH+'_R.txt')
	M_dict['A'] = np.genfromtxt(PATH+'_A.txt')
	return M_dict
	
def ExperimentParamStringSplit(Path):
	#data\exp\INTENSITY.WAVELENGTH.PULSE.ANGLE.SAMPLE.txt
	dummy = Path.split("/")
	List = dummy[2].split(".")
	
	Intensity = int(List[0])
	Wavelength = int(List[1])
	Pulse = int(List[2])
	Angle = int(List[3])
	Sample = int(List[4])
	
	return Intensity, Wavelength, Pulse, Angle, Sample
	
def ExperimentParamPrinter(EXPERIMENT_PATH, I, Parameters, CHECK=1):
	Experiment_Parameters = ExperimentParamStringSplit(EXPERIMENT_PATH)
	print("Experimental Intensity:", Experiment_Parameters[0], "Gw/cm\u00b2")
	print("Model Intensity:", I, "Gw/cm\u00b2 \n")
	
	if Experiment_Parameters[0] != I and CHECK == 1:
		raise ValueError("Experimental intensity does not match with model intensity!")

	print("Experimental Wavelength:", Experiment_Parameters[1], "nm")
	print("Model Wavelength:", int(Parameters['wavelength']/nm), "nm\n")

	if Experiment_Parameters[1] != Parameters['wavelength']/nm and CHECK == 1:
		raise ValueError("Experimental wavelength does not match with model wavelength!", int(Parameters['wavelength']), Experiment_Parameters[1])

	print("Experimental Pulse:", Experiment_Parameters[2], "fs")
	print("Model Pulse:", Parameters['pulse']*1e15, "fs\n")
	
	if Experiment_Parameters[2] != Parameters['pulse']*1e15 and CHECK == 1:
		raise ValueError("Experimental pulse time does not match with model pulse time!")
		
	print("Experimental Angle:", Experiment_Parameters[3], "degrees")
	print("Model Angle:", Parameters['angle'], "degrees\n")

	if Experiment_Parameters[3] != Parameters['angle'] and CHECK == 1:
		raise ValueError("Experimental angle does not match with model angle!")

	print("Experimental Sample Thickness:", Experiment_Parameters[4], "nm")
	print("Model Sample Thickness:", Parameters['thick'][1], "nm")

	if Experiment_Parameters[4] != Parameters['thick'][1] and CHECK == 1:
		raise ValueError("Experimental sample does not match with model sample!")
	
	if CHECK == 1:
		print("\nExperimental parameters equal model inputs.")
	if CHECK == 0:
		print("\nExperimental parameters and model input have not been checked. Please perform this manually.")
	return
