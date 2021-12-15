import sys 
import os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.plotfiles as da

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

da.import_file_as_module( os.path.join(inputFilePath, inputfileName), 'inputfile')
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
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "font.size": 9,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7
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

#----------------------------------------
# start to make figure 
f, axs = plt.subplots(nrows=2, ncols=1, constrained_layout=True)# figsize=(23/2,20/2) )
f.suptitle('Plenum Stats')

# what we want to plot
plotKeyList = ['Pressure', 'Temperature']
# corresponding data units
data_units = [' [bar]', ' [K]']

meanValueMinTime = 300.

for i, ax, subKey, unit in zip(range(len(plotKeyList)), axs.flatten(), plotKeyList, data_units):
   statsData, headerDict = readPlenumStats(output_path, subKey)
   
   #scale Temperature to physical values
   if 'Temperature' in subKey:
      mask = np.ones(statsData.shape, bool)
      mask[:, headerDict['time']] = False # except the time values
      statsData[mask] *= T_ref 

   time = statsData[:, headerDict['time']]
   min = statsData[:, headerDict['min']]
   max = statsData[:, headerDict['max']]
   mean = statsData[:, headerDict['mean']]

   ax.plot(time, mean, color=colors[i], label=subKey)
   ax.fill_between(time, max, min, alpha=0.5, color=colors[i])
   
   meanValue = np.mean(mean[time>meanValueMinTime])
   ax.set_title('mean value = {} {} for '.format(round(meanValue, 2), unit)
               + r'$t \in [{},{}]$'.format(np.rint(meanValueMinTime), np.rint(time[-1]))
               )
   ax.set(xlabel='time', ylabel=subKey+unit, xlim=(time[0], None))
   ax.grid(True)

f.savefig('{}/Plenum_stats.png'.format(output_path), bbox_inches='tight')

f.clear()
plt.close(f)
