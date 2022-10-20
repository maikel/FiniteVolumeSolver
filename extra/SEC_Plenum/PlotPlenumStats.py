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
               +'\t1. argument must be dataPath!\n'
               +'\toptional argument is name of the inputfile\n'
               +'\te.g. {} path --config=inputfile.py'.format(sys.argv[0])
            )
  raise RuntimeError(errMsg)

# parsing the datapath from terminal
dataPath = str(sys.argv[1]) # path to data
if not os.path.exists(dataPath):
   raise FileNotFoundError('given Path: {} does not exist!'.format(dataPath))
inputFilePath = dataPath # assumes inputfile is located in datapath

# name of the inputfile is optional
optional = [ int(el.rsplit('=',1)[-1]) for el in sys.argv if '--config=' in el ]
if not optional:
    optional = ['inputfile.py'] # default value 
inputfileName = optional[0]

import_file_as_module( os.path.join(inputFilePath, inputfileName), 'inputfile')
from inputfile import t_ref, T_ref, ControlOptions

from scipy.signal import savgol_filter
window_length=51 # length of the filter window 51
polyorder=3 # order of the polynomial 
SMOOTHDATA=False

# ftypes = ['.png', '.pdf']
ftypes = ['.png']#, '.pdf']
#--------------------------------
# style settings for matplotlib
plt.style.use('seaborn')
tex_fonts = {
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

usetex = matplotlib.checkdep_usetex(True)
if usetex:
  # Use LaTeX to write all text
  tex_fonts.update({"text.usetex": True, 
                  "font.family": "serif"}) 
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

def readPlenumStatsUnorderd(path, string):
   fname = os.path.join(path, 'Plenum_{}_stats_unorderd.dat'.format(string))
   if not os.path.isfile(fname):
      raise FileNotFoundError()
   headerDict = {}
   with open(fname) as f:
      line = f.readline()
      line = line.split('#')[-1]
      headerDict = {key: value for value, key in enumerate(line.split()) }
   dat = np.loadtxt(fname, skiprows=1)
   ind = np.argsort( dat[:,0] ) # sort times col
   dat = dat[ind] # sort with this indices data
   return dat, headerDict

def smooth(data):
   return savgol_filter(data,     
                  window_length=window_length, 
                  polyorder=polyorder, 
                  mode='mirror')

#----------------------------------------
# start to make figure 

# what we want to plot
plotKeyList = ['Pressure', 'Temperature']
axLabelList = ['p', 'T']
data_units = [' [bar]', ' [K]'] # corresponding data units

meanValueMinTime = 300.

# optional slicing in time-dimension
tplotmin = 0.
tplotmax = 500.0

f, axs = plt.subplots(nrows=2, ncols=1, constrained_layout=True)# figsize=(23/2,20/2) )
f.suptitle('Turbine Plenum Stats')
print('Plot MinMax Figure')

for i, ax, subKey, unit, symbol in zip(range(len(plotKeyList)), axs.flatten(), plotKeyList, data_units, axLabelList):
   try: 
      statsData, headerDict = readPlenumStats(output_path, subKey)
   except:
      statsData, headerDict = readPlenumStatsUnorderd(output_path, subKey)
   
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

   # define labels and strings   
   meanValue = np.mean(mean[time>meanValueMinTime])
   meanString = 'mean {}'.format(symbol)
   fillBetweenString = '+-MinMax({})'.format(symbol)
   xlabel = 'time'
   ylabel = subKey+unit
   title = '{} {} {} {} for '.format(meanString, 'approx', round(meanValue, 2), unit.lstrip())
   title+='t in [{}{}]'.format(np.rint(meanValueMinTime), np.rint(time[-1]))
   if usetex:
      meanString = r'$\overline{{{}}}$'.format(symbol)
      fillBetweenString = r'$\pm min({}) \backslash max({})$'.format(symbol,symbol)
      xlabel = r'$t$'+' [-]'
      ylabel = r'${}$'.format(symbol)+unit
      title = '{} {} {} {} for '.format(meanString, r'$\approx$', round(meanValue, 2), unit[2:-1])
      title+= r'$t \in [{},{}]$'.format(np.rint(meanValueMinTime), np.rint(time[-1]))

   ax.plot( time[t_index_array], 
            mean[t_index_array], 
            color=colors[i], 
            label=meanString
            )
   ax.fill_between( time[t_index_array],
            max[t_index_array], 
            min[t_index_array], 
            alpha=0.5, 
            color=colors[i],
            label=meanString+fillBetweenString 
            )

   ax.set_title(title)
   ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(tplotmin, tplotmax))
   ax.grid(True)
   ax.legend(loc='upper left')

   # ax.label_outer()

for ftype in ftypes:
   f.savefig('{}/Plenum_stats-MinMaxValue{}'.format(output_path, ftype), bbox_inches='tight', dpi=600)

f.clear()
plt.close(f)

#-----------------------------------------------------------
# start to make figure 
f, axs = plt.subplots(nrows=2, ncols=1, constrained_layout=True, figsize=(8.0, 5.))
f.suptitle('Turbine Plenum Stats')
skip=5
stdFactor = 1.0

print('Plot std Figure')
for i, ax, subKey, unit, symbol in zip(range(len(plotKeyList)), axs.flatten(), plotKeyList, data_units, axLabelList):
   try: 
      statsData, headerDict = readPlenumStats(output_path, subKey)
   except:
      statsData, headerDict = readPlenumStatsUnorderd(output_path, subKey)
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

   # define labels and strings   
   meanValue = np.mean(mean[time>meanValueMinTime])
   meanString = 'mean {}'.format(symbol)
   fillBetweenString = '+-{}*std({})'.format(int(stdFactor), symbol)
   xlabel = 'time'
   ylabel = subKey+unit
   title = '{} {} {} {} for '.format(meanString, 'approx', round(meanValue, 2), unit.lstrip())
   title+='t in [{}{}]'.format(np.rint(meanValueMinTime), np.rint(time[-1]))
   if usetex:
      meanString = r'$\overline{{{}}}$'.format(symbol)
      fillBetweenString = r'$\pm {} \cdot \sigma({})$'.format(int(stdFactor), symbol)
      xlabel = r'$t$'+' [-]'
      ylabel = r'${}$'.format(symbol)+unit
      title = '{} {} {} {} for '.format(meanString, r'$\approx$', round(meanValue, 2), unit[2:-1])
      title+= r'$t \in [{},{}]$'.format(np.rint(meanValueMinTime), np.rint(time[-1]))
   
   if SMOOTHDATA:
      MeanData = smooth(mean[t_index_array])
   else:
      MeanData = mean[t_index_array]
   
   ax.plot( time[t_index_array],
            MeanData, 
            color=colors[i], 
            label=meanString )
   
   ax.fill_between(time[t_index_array], 
            MeanData+stdFactor*std[t_index_array], 
            MeanData-stdFactor*std[t_index_array],
            alpha=0.5, 
            color=colors[i],
            label=meanString+fillBetweenString )
   
   meanValue = np.mean(mean[time>meanValueMinTime])
   ax.set_title(title)
   ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(tplotmin, tplotmax)         )
   ax.grid(True)
   ax.legend(loc='upper left')
   # if 'Temperature' in subKey:
   #    ax.set( ylim=(2000,4000) )
   
   # ax.label_outer()

for ftype in ftypes:
   f.savefig('{}/Plenum_stats-std{}'.format(output_path, ftype), bbox_inches='tight', dpi=600)

f.clear()
plt.close(f)