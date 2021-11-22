import matplotlib.pyplot as plt

from IDRIM.constants import *
from IDRIM.relations import *
from IDRIM.commons import *

def M_GraphStaple():
	'''Sets settings for all graphs in notebook. No args. Edit graph.py to change.'''
	textsize = 15
	font = {'family': 'DejaVu Sans', 'weight': 'normal','size': textsize,}
	plt.rc('font', family='DejaVu Sans', serif='Times')
	plt.rc('text', usetex=False)
	plt.rc('xtick', labelsize=textsize)
	plt.rc('ytick', labelsize=textsize)
	plt.rc('axes', labelsize=textsize)
	plt.rc('lines', linewidth=2)
	plt.rcParams["figure.figsize"] = (8, 4)
	width = 10
	height = width / 1.618
	return


def GraphDispersion(TCO, Parabolic):
	'''Plot dispersion relations of TCO and Parabolic. Requires arrays of TCO dispersion and Parabolic.'''
	plt.plot(TCO, energy_array/e, label="TCO")
	plt.plot(Parabolic, energy_array/e, label="Parabolic")
	plt.xlabel("Wave Vector ($m^{-1}$)")
	plt.ylabel("Energy (eV)")
	plt.legend()
	return
	
def GraphDOS(TCO, Parabolic):
	'''Plot density of state relations of TCO and parabolic. Requires arrays.'''
	plt.plot(TCO*e, energy_array/e, label="TCO")
	plt.plot(Parabolic*e, energy_array/e, label="Parabolic")
	plt.xlabel("DOS")
	plt.ylabel("Energy (eV)")
	plt.legend()
	return
	
def GraphMu(Mu):
	'''Plot chemical potential against temperature. Requires an array.'''
	plt.plot(temperature_array, Mu/e)
	plt.xlabel("Temperature (K)")
	plt.ylabel("Chemical Potential (eV)")
	plt.ylim(-4, 2)
	plt.xlim(0, temp_max)
	return
	
def GraphPhononHC(Cp, Xlimit=1200):
	'''Plots phonon heat capacity against temperature. Requires an array. Optional arg sets the limit of x axis. for full graph, set to 20000'''
	plt.plot(temperature_array, Cp)
	plt.xlabel("Temperature (K)")
	plt.ylabel("Phonon Heat Capacity (J/K)")
	plt.xlim(300, Xlimit) #Plataeus after default of 1200
	return

def GraphElectronHC(Ce):
	'''Plot electron heat capacity. Requires an array.'''
	plt.plot(temperature_array[:-1], Ce) #temperature_array slice is needed as differentiation reduces array size by 1
	plt.xlabel("Temperature (K)")
	plt.ylabel("Electron Heat Capacity (J/K)")
	return

def GraphEffectiveMass(Ext, Opt, Con, Wei):
	'''Plots effective mass definitions against energy. Requires arrays.'''
	plt.plot(Con/m_e, energy_array/e, label="Conventional")
	plt.plot(Opt/m_e, energy_array/e, label="Optical", color="green")
	plt.plot(Ext/m_e, energy_array/e, label="Extended", linestyle="dashed", color="black")
	plt.plot(Wei/m_e, energy_array/e, label="Weighted")
	plt.legend()
	plt.xlabel("Effective Mass (/$m_e$)")
	plt.ylabel("Energy (eV)")
	plt.xlim(0,5)
	return

def GraphAvgEffMass(AEM):
	'''Plots average effective mass against temperature. Requires array.'''
	plt.plot(temperature_array, AEM/m_e)
	plt.ylabel("Average Effective Mass (/$m_e$)")
	plt.xlabel("Temperature (K)")
	plt.xlim(0, temp_max)
	return

def GraphPlasmaFrequency(wp):
	'''Plots plasma frequency. Requires array.'''
	plt.plot(temperature_array, wp)
	plt.ylabel("Plasma Frequency")
	plt.xlabel("Temperature (K)")
	return

def GraphRefractiveIndex(RI):
	'''Plots real and imaginary refractive index against temperature. Requires just 1 complex array.'''
	plt.plot(temperature_array, RI.real, label="Real Index")
	plt.plot(temperature_array, RI.imag, label="Imaginary Index")
	plt.legend()
	plt.xlabel("Temperature (K)")
	plt.ylabel("Refractive Index")
	return
	
def Graphing_TimeVSTemperature(Electron_Temp, Phonon_Temp):
	'''Plots time vs temperature plot from model. Requires Electron and Phonon temps from model.'''
	plt.plot(time_array/1e-12, Electron_Temp, label="Electron temperature")
	plt.plot(time_array/1e-12, Phonon_Temp, label="Phonon temperature")
	plt.xlabel("Time (ps)")
	plt.ylabel("Temperature (K)")
	plt.legend()
	return

def Graphing_IntensityVSMaxTemp(MaxTemp):
	'''Plots intensity vs maximum temperature from model. Requires max electron temp from model.'''
	plt.plot(Intensity_array, MaxTemp)
	plt.xlabel("Intensity (GW/cm$^2$)")
	plt.ylabel("Max Temperature (K)")
	return

def Graphing_TemperatureVSCoefficient(Electron_Temp, Transmission, Reflectance, Absorbance):
	'''Plots temperature vs optical coefficients. Requires arrays from model.'''
	plt.plot(Electron_Temp, Transmission, label="Transmission", color="blue")
	plt.plot(Electron_Temp, Reflectance, label="Reflection", color="orange")
	plt.plot(Electron_Temp, Absorbance, label="Absorbance", color="green")
	plt.xlabel("Temperature (K)")
	plt.ylabel("Optical Coefficient")
	plt.legend()
	return

def Graphing_TimeVSCoefficients(Electron_Temp, Transmission, Reflectance, Absorbance):
	'''Plots time vs optical coefficients. Requires arrays then performs an interpolation.'''
	
	plt.plot(time_array/1e-12, Transmission, label="Transmission", color="blue")
	plt.plot(time_array/1e-12, Reflectance, label="Reflectance", color="orange")
	plt.plot(time_array/1e-12, Absorbance, label="Absorbance", color="green")
	plt.xlabel("Time (ps)")
	plt.ylabel("Optical Coefficient")
	plt.legend()
	return

def Graphing_IntensityVSRefractiveIndex(TeMax, Temperature, Refractive_Index):
	'''Plots intensity vs refractive index. Requires Max electron temperature and electron temperature and refractive index. Uses interpolation.'''
	IDRI_array = Interpolate(TeMax, Temperature, Refractive_Index)
	plt.plot(Intensity_array, IDRI_array.real, label="Real")
	plt.plot(Intensity_array, IDRI_array.imag, label="Imaginary")
	plt.xlabel("Intensity GW/cm$^2$")
	plt.ylabel("Refractive index, n")
	plt.legend()
	return

def Graphing_IntensityVSCoefficient(TeMax, Temperature, Refractive_Index):
	'''Plots intensity vs optical coefficients. Requires max electron temp, electron temp and refractive index.'''
	Transmission_Intensity_array = ([])
	Reflectance_Intensity_array = ([])
	Absorbance_Intensity_array = ([])
	
	IDRI_array = Interpolate(TeMax, Temperature, Refractive_Index)
	I_points = int((Intensity_max/Intensity_resolution)+1)
	for i in range(I_points): #Intensity loop begin
		n = IDRI_array[i]
		R, T, A = TMM_Run(n)
		Transmission_Intensity_array.append(T)
		Reflectance_Intensity_array.append(R)
		Absorbance_Intensity_array.append(A)
	#Intensity loop end
	plt.plot(Intensity_array, Transmission_Intensity_array, label="Transmission")
	plt.plot(Intensity_array, Reflectance_Intensity_array, label="Reflectance")
	plt.plot(Intensity_array, Absorbance_Intensity_array, label="Absorbance")
	plt.xlabel("Intensity GW/cm$^2$")
	plt.ylabel("Coefficients")
	plt.legend()
	return
