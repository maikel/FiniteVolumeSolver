import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('extra')[0]
sys.path.append(FVS_path+'/extra/')

import amrex.h5_io as io
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

os.makedirs(output_path, exist_ok=True)

filename_basic = "{}/ControlState.h5".format(dataPath)
print("Read in data from {}".format(filename_basic))

if RESTARTEDSIMULATION:
   datas, times, datas_dict = io.h5_load_restartedTimeseries(filename_basic)
else:
   datas, times, datas_dict = io.h5_load_timeseries(filename_basic)


def getSubKeyList(substring):
   return [ key for key in datas_dict.keys() if substring in key ]

def getCarnotEff(comp_temp, normalized=True, T_ref=T_ref):
   if normalized:
      carnot_eff = 1.0-1.0/comp_temp
   else:
      carnot_eff = 1.0-T_ref/comp_temp
   return carnot_eff

def getSECModeInitTimeID(times, comp_pressure, target_pressure, abs_time_offset=5., relTol=1.0e-4):
   # compute relative pressure difference
   rel_diff = (target_pressure - comp_pressure)/target_pressure
   # get first time where rel_diff is lower than given tolerance
   SECInitId = np.argmax(rel_diff<relTol)
   SECInitTime = times[SECInitId]+abs_time_offset
   print("SEC Mode was initiated at time = {}".format(SECInitTime-abs_time_offset))
   print("For Post Processing absolute time offset = {} was used!".format(abs_time_offset))
   SECInitIds = times>SECInitTime
   return SECInitIds, SECInitTime
#----------------------------------------
# start to make figure 
f, axs = plt.subplots(nrows=4, ncols=2, figsize=(23/2,20/2) )
f.suptitle('Control Station')

# what we want to plot
plotKeyList = ['compressor_mass_flow', 'turbine_mass_flow', 'power', 'rpm', 'pressure', 'temperature', 'fuel_consumption_rate', 'efficiency']
# corresponding data units
data_units = ['', '', '', ' [-]', ' [bar]', ' [K]', '', '', '']

# normalize rpm with max value
datas[:, datas_dict['current_rpm']] = datas[:, datas_dict['current_rpm']] / ControlOptions['rpmmax']
print( "maximum compressor_pressure = {}".format(np.max(datas[:, datas_dict['compressor_pressure']])) )

# list of keys where PLOTFILLBETWEEN is applied:
plotbetw_list = ['compressor_mass_flow_in', 'compressor_mass_flow_out', 'turbine_mass_flow_in', 'turbine_mass_flow_out', 'compressor_power', 'power_out', 'turbine_power', 'turbine_pressure', 'turbine_temperature', 'efficiency', 'fuel_consumption_rate']

t_id_sec_init, t_sec_init = getSECModeInitTimeID(times, datas[:, datas_dict['compressor_pressure']], ControlOptions['target_pressure_compressor'])

if t_sec_init>times[-1]:
   PLOTFILLBETWEEN=False
   print("PLOTFILLBETWEEN was set to False, because t_sec_init={} < tend={}".format(t_sec_init, times[-1]))

if PLOTFILLBETWEEN:
   FVS_skip=1
else:
   FVS_skip = 389 # skip time data to make plot clearer # 389
datas = datas[::FVS_skip, :]
times = times[::FVS_skip]
print("FVS delta t is {}".format(np.diff(times)[0]))

for i, ax, subKey in zip(range(len(plotKeyList)), axs.flatten(), plotKeyList):
   subKeyList = getSubKeyList(subKey)
   for j, key in enumerate(subKeyList):
      if i<2 or i>6: # don't destroy the labels from the last two plots
         lab = key.replace('_', ' ')
      else:
         lab = key.replace(subKey, '')
         lab = lab.replace('_', ' ')

      if (i<2):
         if ('comp' in lab):
            lab = lab.replace('compressor', 'comp. plenum')
            lab = lab.replace('mass flow', '')
         elif ('turb' in lab):
            lab = lab.replace('turbine', 'turb. plenum')
            lab = lab.replace('mass flow', '')
      
      if 'temperature' in key:
         datas[:, datas_dict[key]] *=T_ref
      
      if USEONLYSAVGOLFILTER:
         if ('SEC_Mode' in subKey) or ('rpm' in subKey):
            ax.plot( times, datas[:, datas_dict[key]], label=lab )
         else:
            cutoff=-1 #-window_length//2
            ax.plot( times, savgol_filter(datas[:, datas_dict[key]], window_length=window_length, polyorder=polyorder, mode='mirror'), 
                        color=colors[j], label=lab )
            # ax.plot( times[cutoff:], datas[cutoff:, datas_dict[key]] )
      elif PLOTFILLBETWEEN:
         if key in plotbetw_list:
            # first plot solid lines before SEC was stable
            ax.plot( times[~t_id_sec_init], datas[~t_id_sec_init, datas_dict[key]], color=colors[j], label=lab )
            # after SEC is stable plot mean value and shaded areas
            dat = datas[:, datas_dict[key]]
            peaks, _ = find_peaks(dat) # get local maxima, peaks are the indices from the local maxima in the 1d array
            low_peaks, _ = find_peaks(-dat) # get local minima (or maxima from mirrored data on the x axis)

            yhi = np.interp(times[t_id_sec_init], times[peaks], dat[peaks])
            ylo = np.interp(times[t_id_sec_init], times[low_peaks], dat[low_peaks])

            # ax.plot( times[t_id_sec_init], np.mean([yhi,ylo], axis=0), color=colors[j]) # plot mean value for nonoscillating data
            
            # plot mean value for oscillating data
            skip=389
            skip_t = times[::skip]
            skip_dat = dat[::skip]
            ax.plot( skip_t[t_id_sec_init[::skip]], savgol_filter(skip_dat[t_id_sec_init[::skip]], window_length=window_length, 
            polyorder=polyorder, mode='nearest'), color=colors[j]) 

            # plot shaded area
            ax.fill_between(times[t_id_sec_init], yhi, ylo, alpha=0.5, color=colors[j]) 
         else:   
            ax.plot( times, datas[:, datas_dict[key]], color=colors[j], label=lab )
      else:
         ax.plot( times, datas[:, datas_dict[key]], color=colors[j], label=lab )

      if i==3 or i>=6:
         pass
      else:
         ax.legend(loc='upper left')
      if 'efficiency' in subKey:
         ax.plot(times, getCarnotEff(datas[:, datas_dict['compressor_temperature']], normalized=False), color=colors[1], label='Carnot efficiency')
         ax.legend(loc='lower right')
         ax.set(ylim=(0,0.5))
      
      if 'flow' in subKey:
         ylab = subKey.replace('_', ' ')
         ylab = ylab.split(' ', 1)[-1]
      else:
         ylab = subKey.replace('_', ' ')

      if 'rpm' in ylab:
         ylab='rpm / rpmmax'
      ylab += data_units[i]
         
      ax.set(xlabel='time', ylabel=ylab, xlim=(times[0], None))
      ax.grid(True)

ax = axs.flatten()
# ax[1].set(ylim=(0, None))

# set mass flow plots to same ylim!
comp_ylim = ax[0].get_ylim()
turb_ylim = ax[1].get_ylim()
if comp_ylim[1]<turb_ylim[1]:
   ax[0].set_ylim(turb_ylim)
else:
   ax[1].set_ylim(comp_ylim)

f.subplots_adjust(hspace=0.3, wspace=0.3)
f.savefig('{}/control_state.png'.format(output_path), bbox_inches='tight')
# f.savefig('{}/control_state.pdf'.format(output_path), bbox_inches='tight')

f.clear()
plt.close(f)
