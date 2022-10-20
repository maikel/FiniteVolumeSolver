import math
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import json

#-------------------------------
# dicts
valueDict = {
  'time' : {
    'data' : None,
    'label': 'time',
    'symbol': 't',
    'latexLabel': r'$t$',
    'dim': '[s]'
  },
  'pressure' : {
    'data' : None,
    'label': 'pressure',
    'symbol': 'p',
    'latexLabel': r'$p$',
    'dim': '[bar]'
  },
  'density' : {
    'data' : None,
    'label': 'density',
    'symbol': 'rho',
    'latexLabel': r'$\rho$',
    'dim': r'$[kg/m^3]$'
  },
  'temperature' : {
    'data' : None,
    'label': 'temperature',
    'symbol': 'T',
    'latexLabel': r'$T$',
    'dim': '[K]'
  },
  'specificVolume' : {
    'data' : None,
    'label': 'specific Volume',
    'symbol': 'v',
    'latexLabel': r'$v$',
    'dim': r'$[m^3/kg]$'
  },
  'specificEntropy' : {
    'data' : None,
    'label': 'specific Entropy difference',
    'symbol': 's',
    'latexLabel': r'$\Delta s$',
    'dim': r'$[J/(kg\cdot K)]$'
  },
  'specificEnthalpy' : {
    'data' : None,
    'label': 'specific Enthalpy',
    'symbol': 'h',
    'latexLabel': r'$h$',
    'dim': '[J/kg]'
  },
  'fuel' : {
    'data' : None,
    'label': 'fuel',
    'symbol': 'Yf',
    'latexLabel': r'$Y_F$',
    'dim': '[-]'
  },
}
#----------------------------------------------
## classes

# fucking awesome can hide print statements
# usage see class below
# https://stackoverflow.com/a/45669280
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# reference parameter
class pv_diagramm_class(object):
    def __init__(self):  
        with HiddenPrints():
            # perfect gas constants
            self.specGasConstant = 1.0 # [-]
            self.gamma = 1.4 # [-]

            # reference parameter
            self.ref_length = 1.0 # [m]
            self.ref_pressure = 1.0e5 # [Pa]
            self.ref_temperature = 300. # [K]
            self.ref_GasConstant = 287.4 # [J/kg/K] should be dry air
            self.ref_density = self.ref_pressure / self.ref_temperature / self.ref_GasConstant # [kg/m^3]
            self.ref_velocity = math.sqrt(self.ref_pressure / self.ref_density) # [m/s]
            self.ref_time = self.ref_length / self.ref_velocity # [s]

            # dimless ambientConditions
            self.ambientPressure = 1.0e5 / self.ref_pressure
            self.ambientTemperature = 273.15 / self.ref_temperature
            self.ambientDensity = self.computeDensity(self.ambientPressure, self.ambientTemperature)
    
    def prettyprint(self):
      return json.dumps(vars(self), sort_keys=False, indent=4)

    def __setattr__(self, __name, __value) -> None:
        print('Please use class method setAttribute to change private attributes')
        super().__setattr__(__name, __value)

    # please use this to change private class attributes
    def setAttribute(self, name, value):
        if not hasattr(self, name):
            raise AttributeError("this class has no atttribute named %s !"%name)
        else:
            super().__setattr__(name, value)

    def getSpecGasConstant(self, Dimless=True):
        Rspec = self.specGasConstant
        if not Dimless:
            Rspec *= self.ref_GasConstant
        return Rspec
    
    def getHeatCapacityConstantVolume(self, Dimless=True):
        gammaMinusOne = self.gamma-1.
        Rspec = self.getSpecGasConstant(Dimless)
        return Rspec/gammaMinusOne

    def getHeatCapacityConstantPressure(self, Dimless=True):
        gammaMinusOne = self.gamma-1.
        Rspec = self.getSpecGasConstant(Dimless)
        return Rspec*self.gamma/gammaMinusOne
    
    def computeEntropy(self, T, rho, s0=0., Dimless=True):
        cv = self.getHeatCapacityConstantVolume(Dimless)
        Rspec = self.getSpecGasConstant(Dimless)
        s = np.zeros_like(rho)
        s[0] = s0 # set first element
        sdiff = (cv * np.log(T[1:]/T[:-1]) 
                    + Rspec * np.log(rho[:-1]/rho[1:])
                )
        for i in range(s.shape[0]-1):
            s[i+1] = s[i]+sdiff[i]
        return s

    def computeEnthalpy(self, T, Dimless=True):
        # H/m = h = cp * T
        cp = self.getHeatCapacityConstantPressure(Dimless)
        return cp*T

    #-----------------------------------------------------------------
    def computeTemperature(self, p_dimless, rho_dimless, Dimless=True):
        # p_dimless and rho_dimless are assumed to be Dimless!
        # in self.ref_temperature is self.ref_GasConstant included!!
        ## --> so we need Rspec always Dimless
        Rspec = self.getSpecGasConstant(True)
        if Dimless:
          return p_dimless / rho_dimless / Rspec
        else:
          return p_dimless / rho_dimless / Rspec * self.ref_temperature

    def computeDensity(self, p_dimless, T_dimless, Dimless=True):
        Rspec = self.getSpecGasConstant(True)
        if Dimless:
          return p_dimless/T_dimless/Rspec
        else: 
          return p_dimless/T_dimless/Rspec*self.ref_density

    def computePressure(self, rho_dimless, T_dimless, Dimless=True):
        Rspec = self.getSpecGasConstant(True)
        if Dimless:
          return (rho_dimless*T_dimless*Rspec)
        else:
          return (rho_dimless*T_dimless*Rspec)*self.ref_pressure
    #--------------------------------------------------------------


#-------------------------------------
# functions

def getNumpyMachineEpsilon():
  return np.finfo(float).eps

def findNearest1D(array, value):
  array = np.asarray(array)
  sortedIndices = np.argsort(np.abs(array - value))
  
  idPos = []
  idNeg = []
  # greater than value
  for id in sortedIndices:
    if array[id]>=value:
      idPos.append(id)
      break
  # less than value
  for id in sortedIndices:
    if array[id]<value:
      idNeg.append(id)
      break
  if (not idPos):
    raise AssertionError(f"no value is greater than {value}!")
  elif (not idNeg):
    raise AssertionError(f"no value is less than {value}!")
  
  # get rid of list type
  idPos, = idPos
  idNeg, = idNeg

  return idNeg, idPos

def findNearest2D(array, value):
  absArray = np.abs(array-value)
  index = np.where( absArray == np.amin(absArray) )
  return index

def checkMin(value, array):
  amin = np.amin(array)
  if (value < amin):
    # print('scalar value: {} is lower than min scalar value: {}'.format(value, amin))
    return False
  else:
    return True

def checkMax(value, array):
  amax = np.amax(array)
  if (value > amax):
    # print('scalar value: {} is greater than max scalar value: {}'.format(value, amax))
    return False
  else:
    return True


def valueCheck(value, array):
    MIN = checkMin(value, array)
    MAX = checkMax(value, array)
    if not MIN:
      return False
    if not MAX:
      return False
    return True

def interpolate1D(value, x, y, ids):
      data = np.interp( value, 
                [x[ids[0]], x[ids[1]]], 
                [y[ids[0]], y[ids[1]]] 
               )
      return data

def convertStringToFloat(string):
  try:
    return float(string)
  except ValueError:
    return string

def getRelevantControlStateData(timepoint, csData, csTimes):
  return csData[csTimes>=timepoint, :][0]
#----------------------------------
# plot routines

def quiverPlot(strX, strY, valueDict, ax, location=[], sl=(), usetex=False):
  # print(f'Plotting {strY}-{strX} diagram')
  color = ['m', 'k', 'g', 'tab:orange', 'c']
  kwargs = {
    'scale_units' : 'xy', 
    'angles' : 'xy', 
    'scale' : 1,
    'zorder' : 7,
  }
  if location:
    # print(f'DEBUG: locations we have seen are {set(location)}')
    kwargs['color'] = [color[i] for i in location]
  
  x = valueDict[strX]['data']
  y = valueDict[strY]['data']

  ax.set(
    xlim = (math.floor(1.1*x.min()), round(1.1*x.max(), 1)),
    ylim = (math.floor(1.1*y.min()), round(1.1*y.max()))
  )

  if sl:
    x = x[sl]
    y = y[sl]
  
  ax.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1], 
                **kwargs)
  key = 'label'
  if usetex:
    key = 'latexLabel'
  ax.set(
      xlabel=f"{valueDict[strX][key]} {valueDict[strX]['dim']}", 
      ylabel=f"{valueDict[strY][key]} {valueDict[strY]['dim']}"
      )

def simplePlot(strX, strY, path, test_scalar, markerstyle='x', tube_id=0, usetex=False):
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(valueDict[strX]['data'], valueDict[strY]['data'], 
              '-', marker=markerstyle)
    key = 'label'
    if usetex:
      key = 'latexLabel'
    ax.set(
      xlabel=f"{valueDict[strX][key]} {valueDict[strX]['dim']}", 
      ylabel=f"{valueDict[strY][key]} {valueDict[strY]['dim']}"
      )
    fig.savefig('{}/{}-{}{}Diagramm_tubeID{}.png'.format(path, str(int(test_scalar)).zfill(5), valueDict[strY]['symbol'], valueDict[strX]['symbol'], tube_id), bbox_inches='tight')

def computeIsobarCurveTemperature(T1, v1, specVolArray):
  # return TemperatureArray
  const = v1/T1
  return const / specVolArray

def computeIsobarCurveSpecVol(T1, v1, temperatureArray):
  # return TemperatureArray
  const = v1/T1
  return const * temperatureArray

def computeIsochorCurvePressure(p1, T1, temperatureArray):
  # return pressureArray
  const = p1/T1
  return const * temperatureArray

def computeIsochorCurveTemperature(p1, T1, pressureArray):
  # return temperatureArray
  const = p1/T1
  return pressureArray / const

def computeIsothermCurve(p1, v1, isoVolArray):
  eps = getNumpyMachineEpsilon() # prevent zerodivision warning
  const = p1 * v1 # p * v = const
  return const/(isoVolArray+eps)

def computeIsentropCurvePressure(p1, v1, gamma, isoVolArray):
  eps = getNumpyMachineEpsilon() # prevent zerodivision warning
  const = p1 * v1**gamma # p * v^gamma = const
  return const/((isoVolArray+eps)**gamma)

def computeIsentropCurveIsoVol(p1, v1, gamma, pressureArray):
  eps = getNumpyMachineEpsilon() # prevent zerodivision warning
  const = p1 * v1**gamma # p * v^gamma = const
  return (const/(pressureArray+eps))**(1./gamma)

def computeHugoniotCurve(p1, v1, gamma, pressureArray):
  eps = getNumpyMachineEpsilon() # prevent zerodivision warning
  gammaMinusOne = gamma-1.
  pterm = 0.5*(pressureArray+p1+eps)
  counter = p1/gammaMinusOne + pterm
  denominator = pressureArray/gammaMinusOne + pterm
  return v1 * counter / denominator

def drawIsothermLine(pressurePoint, specVolPoint, ax):
  isoVol = np.linspace(*ax.get_xlim(), 101 )
  pressure = computeIsothermCurve(pressurePoint, specVolPoint, isoVol)
  ax.plot(isoVol, pressure, 'r--', label='isotherm '+r'$p\cdot v$', zorder=5)

def drawIsentropLine(pressurePoint, specVolPoint, gamma, ax):
  isoVol = np.linspace(*ax.get_xlim(), 101 )
  pressure = computeIsentropCurvePressure(pressurePoint, specVolPoint, gamma, isoVol)
  ax.plot(isoVol, pressure, 'b--', label='isentrop '+r'$p\cdot v^\gamma$', zorder=5)

def drawHugoniotCurve(pressurePoint, specVolPoint, gamma, ax):
  pressureArray = np.linspace(*ax.get_ylim(), 101 )
  isoVol = computeHugoniotCurve(pressurePoint, specVolPoint, gamma, pressureArray)
  ax.plot(isoVol, pressureArray, 'g--', label='hugoniot curve', zorder=5)

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), loc='best')

# find the point when pressure is rising
def findPressureRisingPoint(pressure, specVol):
  index = np.argmax(np.diff(pressure)>0)
  return (pressure[index], specVol[index])