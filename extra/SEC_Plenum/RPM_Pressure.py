import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import math


R_ref = 287.
p_ref = 101325.
T_ref = 300.
L_ref = 1.0
rho_ref = p_ref / T_ref / R_ref
u_ref = math.sqrt(p_ref / rho_ref)
t_ref = L_ref / u_ref


rpmmin = 10000.0 / 60.0 * t_ref
rpmmax = 60000.0 / 60.0 * t_ref
rpmmean = 0.5 * (rpmmax + rpmmin)
pratiomin = 1.0
pratiomax = 50.0
pratiomean = 0.5 * (pratiomax + pratiomin)
pratiovar = pratiomax - pratiomin
c_0 = 22.0 / 17.0

def CompressorPressureRatio(rpm):
   pressure = pratiomean + c_0 * pratiovar * np.arctan( np.pi * (rpm / rpmmean - 1.0) ) / np.pi
   return pressure

rpm = np.linspace(rpmmin,rpmmax, 100)

f, ax = plt.subplots()

ax.plot(rpm/rpmmax, CompressorPressureRatio(rpm))
ax.set(xlabel='rpm / rpmmax', ylabel='p / p_ref', xlim=(0.1, 0.4), ylim=(0,10))
ax.grid(True)
f.savefig('RPM_pressure_correlation.png', bbox_inches='tight')
