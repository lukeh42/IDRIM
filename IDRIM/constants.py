import numpy as np #Needed for array generation
from scipy.constants import * #Needed for some definitions

#Physical Constants
Angstrom = 1e-10
nm = 1e-9
T0 = 300 #Ambient Temperature (K)
Gwcm2 = 1e+13 #Power Units

#Material Constants !M! Note that these can and will change with ITO samples! 
#These constants were calculated for a 407nm sample.
N =  9.531746606764642e+26			#Electron Density (per m^3) old file had 2.0585e+27
eps_inf = 3.45						#High Frequency Permittivity Limit
me_min = (0.2623384490557345*m_e)	#Effective Mass at Bottom of Band (old 0.3964)
eta = (0.4191/e) 					#Non-parabolic constant. Old file was C #changed for clarity
wp0 = (2.5325e15)		 				#Plasma Frequency at rest.
gamma0 = 291e12						#Drude Scattering Rate at 300K (Assumed to be a constant)

#These constants shouldn't change
TD = 1000 							#Debye Temperature (K)
Ef = 0.96*eV							#Fermi Level at 300K was 1.036327*eV #Now is the Fermi Energy
vf = np.sqrt(2*Ef/me_min) 			#Fermi Velocity at 300K
Tf = Ef/k							#Fermi Temperature (at 1eV, this is 11,600K)
Lf = 1.577483e-14*vf 				#Mean Free Path (Relaxation rate at 300K * vf)
n_substrate = 1.5					#Refractive Index of Substrate

#Convenience constants (Absorptive loss and band non-parabolicity as a physical origin of large nonlinearity in epsilon-near-zero materials. R. Secendo et al.) #not used in main code, still used in some side modules
a = 1/(2*eta)
b = (me_min)/(2*eta*hbar**2)

#System Parameters that shouldn't need to be changed
temp_max = 20000 					# Maximum Temperature (K)
Intensity_max = 400					# Max Intensity, GW/cm^2 #400 as that is the highest anyone has taken sample
DIM = 1000							# Length of arrays
time_points = 50000					# Time Resolution

#System arrays that shouldn't need changing.
time_array = np.linspace(-1500, 1500, time_points)*1e-15	 	# Time array.
energy_array = np.linspace(0, 4*e, DIM)							# Energy array. Setting higher does nothing.
temperature_array = np.linspace(300, temp_max, DIM)				# Temperature array.

Dt = time_array[1]-time_array[0] #saves writing it over and over

Params_path = "data\core\ParameterHistory.txt"
