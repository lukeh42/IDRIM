import numpy as np
from numpy import inf
from scipy.constants import *
import matplotlib.pyplot as plt

from IDRIM.constants import *
from IDRIM.relations import *
from IDRIM.graph import *
from IDRIM.commons import *
from IDRIM.checks import *
from IDRIM.checks import *
from IDRIM.genesis import *
from IDRIM.ymodcom import *

def GlobalSubstrateError(WaveL, Plasma, angle, d_array_end=[10000,14000]):
	'''Inputs: Wavelength, Plasma (non-array), angle.
	Optional: end points of thickness array.
	Returns Delta which is the error.'''
	R = np.zeros(DIM)
	T = np.zeros(DIM)
	A = np.zeros(DIM)
	
	n = np.sqrt(Permittivity(WavelengthToFrequency(WaveL), Plasma))
	
	Thickness_array = np.linspace(d_array_end[0], d_array_end[1], DIM)
	
	for i in range(DIM):
		Thicks = [inf, 407, Thickness_array[i], inf]
		R[i], T[i], A[i] = TMM_Run(n, WaveL, Thicks, angle)
	
	Delta = np.zeros(3)
	#done like this so can have access to all information easily
	Delta[1] = max(T)
	Delta[2] = min(T)
	Delta[0] = Delta[1]-Delta[2] 
	#not divided by 2 as if you are at bottom of peak due to substrate thickness being out of phase, need full delta
	
	for j in range(DIM):
		if T[j] == Delta[1]:#max(T) check
			sub_max  = Thickness_array[j]
		if T[j] == Delta[2]:#min(T) check
			sub_min = Thickness_array[j]
	subs = [sub_max, sub_min]
	
	return Delta, subs
