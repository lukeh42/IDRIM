import numpy as np
from numpy import inf
from scipy.constants import *
import matplotlib.pyplot as plt

from IDRIM.relations import *
from IDRIM.graph import *
from IDRIM.commons import *
from IDRIM.checks import *
from IDRIM.modcom import *
from IDRIM.genesis import *
from IDRIM.ymodcom import *
from IDRIM.file import *
import IDRIM.models.y1 as y


def AllInOne(Parameters, Intensity, MODCON, peak_point=0.0):
	'''
	All In One Model Function.
	Generates needed values from extended drude model, then runs the two temperature model with this.
	Designed to be the end-user function.
	'''
	#Extended Drude Model Steps
	mu_array, Cp_array, Ce_array, wp_array, RI_array, Fit = Regeneration(Parameters, OVERRIDE=1, Text=0)
	Relations_dict = {'mu':mu_array, 'Cp':Cp_array, 'wp':wp_array, 'Ce':Ce_array, 'RI':RI_array, 'Fit':Fit}
	
	#Two-Temperature Model Steps
	Model_Output = y.Core(Intensity, Relations_dict, Parameters, MODCON, peak_point)

	return Model_Output