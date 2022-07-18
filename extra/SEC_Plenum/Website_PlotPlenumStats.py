import sys 
import os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

from amrex.other import import_file_as_module

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# check cli
if len(sys.argv)<2:
   errMsg = ('Not enough input arguments!\n'
               +'\tfirst argument must be dataPath!')
   raise RuntimeError(errMsg)

# parsing the datapath from terminal
dataPath = str(sys.argv[1]) # path to data
if not os.path.exists(dataPath):
   raise FileNotFoundError('given Path: {} does not exist!'.format(dataPath))
inputFilePath = dataPath # assumes inputfile is located in datapath

try:
   inputfileName = str(sys.argv[2]) # optional name of the inputfile
except: 
   inputfileName = 'SEC_Plenum_Arrhenius.py'

import_file_as_module( os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import t_ref, T_ref, ControlOptions

#---------------------------------------
# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = False

from scipy.signal import savgol_filter
USEONLYSAVGOLFILTER=False 
window_length=51 # length of the filter window 51
polyorder=3 # order of the polynomial 

PLOTFILLBETWEEN=True
from scipy.signal import find_peaks
#--------------------------------
# style settings for matplotlib
plt.style.use('seaborn')
tex_fonts = {
    # Use LaTeX to write all text
   #  "text.usetex": True,
   #  "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "font.size": 12,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10
}
plt.rcParams.update(tex_fonts)

# list of named colors which will be used
colors=['b', 'g', 'r', 'm', 'k']
#------------------------------------

outPath = dataPath
output_path = '{}/Visualization'.format(outPath)
print("Processing data in {}".format(outPath))

os.makedirs(output_path, exist_ok=True)

def readPlenumStats(path, string):
   fname = os.path.join(path, 'Plenum_{}_stats.dat'.format(string))
   if not os.path.isfile(fname):
      raise FileNotFoundError()
   headerDict = {}
   with open(fname) as f:
      line = f.readline()
      line = line.split('#')[-1]
      headerDict = {key: value for value, key in enumerate(line.split()) }
   dat = np.loadtxt(fname, skiprows=1)
   return dat, headerDict

def smooth(data):
   return savgol_filter(data,     
                  window_length=window_length, 
                  polyorder=polyorder, 
                  mode='mirror')

#----------------------------------------
# start to make figure 
f, axs = plt.subplots(nrows=2, ncols=1, constrained_layout=True, figsize=(8.0, 5.))# figsize=(23/2,20/2) )
# f.suptitle('Plenum Stats')

# what we want to plot
plotKeyList = ['Pressure', 'Temperature']
axLabelList = ['p', 'T']
# corresponding data units
data_units = [' [bar]', ' [K]']

meanValueMinTime = 300.

# optional slicing in time-dimension
tplotmin = 0.0
# splittedPath = dataPath.rsplit('/',1)[-2] #-2 because / at the end
# if 'vol40.0' in splittedPath:
#    tplotmin = 200.0
# elif 'vol20.0' in splittedPath:
#    tplotmin = 145.0
tplotmax = 500.0
skip=5

for i, ax, subKey, unit, symbol in zip(range(len(plotKeyList)), axs.flatten(), plotKeyList, data_units, axLabelList):
   statsData, headerDict = readPlenumStats(output_path, subKey)
   statsData = statsData[::skip]
   
   #scale Temperature to physical values
   if 'Temperature' in subKey:
      mask = np.ones(statsData.shape, bool)
      mask[:, headerDict['time']] = False # except the time values
      statsData[mask] *= T_ref 

   time = statsData[:, headerDict['time']]
   min = statsData[:, headerDict['min']]
   max = statsData[:, headerDict['max']]
   mean = statsData[:, headerDict['mean']]
   std = statsData[:, headerDict['std']]

   tplotmax = tplotmax if time[-1]>=tplotmax else np.rint(time[-1])
   t_index_array = (time>=tplotmin) & (time<=tplotmax)
   # check if index array is empty
   if not np.any(t_index_array):
      raise IndexError('time index array is empty!')

   meanString = r'$\overline{{{}}}$'.format(symbol)
   ax.plot( time[t_index_array], 
            smooth(mean[t_index_array]), 
            color=colors[i], 
            label=meanString )

   # ax.plot(time[t_index_array], mean[t_index_array], color=colors[i], label=subKey.lower())

   # ax.fill_between(time[t_index_array], max[t_index_array], min[t_index_array], alpha=0.5, color=colors[i])
   
   fac = 1.0
   ax.fill_between(time[t_index_array], 
            smooth(mean[t_index_array])+fac*std[t_index_array], 
            smooth(mean[t_index_array])-fac*std[t_index_array],
            alpha=0.5, 
            color=colors[i],
            label=meanString+r'$\pm\sigma({})$'.format(symbol) )
   
   meanValue = np.mean(mean[time>meanValueMinTime])
   ax.set_title('{} {} {} {} for '.format(meanString, r'$\approx$', int(np.rint(meanValue)), unit[2:-1])
               + r'$t \in [{},{}]$'.format(np.rint(meanValueMinTime), np.rint(time[-1]))
               )
   ax.set(xlabel=r'$t$'+' [-]', 
            ylabel=r'${}$'.format(symbol)+unit, 
            xlim=(tplotmin, tplotmax)
         )
   ax.grid(True)
   ax.legend(loc='best')
   if 'Temperature' in subKey:
      ax.set( ylim=(2000,4000) )
   
   ax.label_outer()

for ftype in ['.png', '.pdf']:
   f.savefig('{}/Plenum_stats-poster{}'.format(output_path, ftype), bbox_inches='tight', dpi=600)

f.clear()
plt.close(f)
