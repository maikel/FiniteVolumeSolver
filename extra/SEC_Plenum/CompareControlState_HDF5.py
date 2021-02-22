import sys, os

# get the absolute path to the FUB FVS-Solver
pathname = os.path.dirname(sys.argv[0])
pathname = os.path.abspath(pathname)
FVS_path = pathname.split('FiniteVolumeSolver')[0]+'FiniteVolumeSolver'
# print(FVS_path) 
sys.path.append(FVS_path+'/extra/')
from amrex.plotfiles import h5_load_timeseries
import amrex.plotfiles as da

import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt


os.environ['HDF5_USE_FILE_LOCKING'] = 'False'

# optional parsing the datapath from the terminal
if (len(sys.argv)>1):
   dataPath = str(sys.argv[1])
   inputFilePath = dataPath
else:
   dataPath = FVS_path+"/build_2D-Release/average_massflow"
   inputFilePath = FVS_path+"/examples/AMReX/EB/2D/"

# bool to read all existing HDF5 files
# this make only sense if we restarted the simulation form the last checkpoint!!
RESTARTEDSIMULATION = True

sys.path.append(inputFilePath)
from SEC_Plenum_Arrhenius import t_ref, ControlOptions, T_ref, rho_ref

outPath = dataPath
output_path = '{}/Visualization'.format(outPath)

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

os.makedirs(output_path, exist_ok=True)

filename_basic = "{}/ControlState.h5".format(dataPath)

if RESTARTEDSIMULATION:
   import glob
   fileNameList = [filename_basic]

   ###### collect data begin
   # check if other h5.* files exist and append them to the list
   otherFiles = glob.glob("{}.*".format(filename_basic))
   if otherFiles:
      fileNameList.append( *otherFiles )

   print(fileNameList)

   # Read in data
   # Attention the last file is the latest!!
   # for example we have Filename.h5 and Filename.h5.1 the last one contains the first data!
   print("Read in data from {}".format(fileNameList[-1]))
   FVS_data, FVS_times, FVS_dict = h5_load_timeseries(fileNameList[-1])

   for filename in reversed(fileNameList[:-1]):
      print("Read in data from {}".format(filename))
      data, time, _ = h5_load_timeseries(fileNameList[0])
      FVS_data = np.concatenate((FVS_data, data))
      FVS_times = np.concatenate((FVS_times, time))
else:
   print("Read in data from {}".format(filename_basic))
   FVS_data, FVS_times, FVS_dict = h5_load_timeseries(filename_basic)

FVS_skip = 400 # skip time data to make plot clearer
FVS_data = FVS_data[::FVS_skip, :]
FVS_times = FVS_times[::FVS_skip]
print("FVS delta t is {}".format(np.diff(FVS_times)[0]))

# get reference data from Kleins Code
GT_filename = "{}/GT_time_series.txt".format(dataPath)

GT_data, GT_times, GT_dict = da.get_controlState_Klein(GT_filename)

GT_skip=14 # skip time data to make plot clearer
GT_data = GT_data[::GT_skip,:]
GT_times = GT_times[::GT_skip]
print("GT delta t is {}".format(np.diff(GT_times)[0]))

# get maximum time from both simulations and create index array
tEnd = min(np.max(FVS_times), np.max(GT_times))
# print(tEnd)
GT_t_index = GT_times<tEnd
FVS_t_index = FVS_times<tEnd
# print(GT_data.shape)
# print(FVS_data.shape)
# print(GT_data[GT_t_index].shape)
# print(FVS_data[FVS_t_index].shape)

################# Calc Mean Values ############################
from scipy.integrate import simps
from scipy.integrate import trapz
def getMeanValue(y, x, mode='scipy_simps'):
   dist = abs(x[-1]-x[0])
   if mode=='numpy_trapz':
      mean = np.trapz(y, x)
      mean *=(1.0/dist)
      return mean
   elif mode=='numpy_mean':
      mean = np.mean(y)
      return mean
   elif mode=='scipy_simps':
      mean = simps(y, x)/dist
      return mean
   elif mode=='scipy_trapz':
      mean = trapz(y, x)/dist
      return mean

def CompareMeanValue(gt_data, gt_time, fvs_data, fvs_time, variable, ndig=4):
   gt_mean = getMeanValue(gt_data, gt_time)
   fvs_mean = getMeanValue(fvs_data, fvs_time)
   means = [ ['GT Mean', 'FVS Mean', 'abs diff'],
         [gt_mean, fvs_mean, abs(gt_mean-fvs_mean)] ]
   means[1] = [ round(val, ndig) for val in means[1] ]

   format_row = '{:>12}'*len(means[0])
   print("Means for {}:".format(variable))
   for row in means:
      print(format_row.format(*row))
   print()
   return means

tInit = 150.0
GT_idx_array = (GT_times>tInit) & (GT_times<tEnd )
FVS_idx_array = (FVS_times>tInit) & (FVS_times<tEnd )

gt_time = GT_times[GT_idx_array]
fvs_time= FVS_times[FVS_idx_array]

meanVariables = []
meanValues = []

var = 'fuel_consumption_rate'
means = CompareMeanValue(GT_data[GT_idx_array, GT_dict[var]], gt_time, FVS_data[FVS_idx_array, FVS_dict[var]], fvs_time, var, ndig=4)
meanVariables.append(var); meanValues.append(means)

var = 'efficiency'
means = CompareMeanValue(GT_data[GT_idx_array, GT_dict[var]], gt_time, FVS_data[FVS_idx_array, FVS_dict[var]], fvs_time, var, ndig=4)
meanVariables.append(var); meanValues.append(means)

var = 'turbine_pressure'
means = CompareMeanValue(GT_data[GT_idx_array, GT_dict[var]], gt_time, FVS_data[FVS_idx_array, FVS_dict[var]], fvs_time, var, ndig=4)
meanVariables.append(var); meanValues.append(means)

var = 'turbine_temperature'
means = CompareMeanValue(GT_data[GT_idx_array, GT_dict[var]], gt_time, FVS_data[FVS_idx_array, FVS_dict[var]], fvs_time, var, ndig=4)
meanVariables.append(var); meanValues.append(means)

meanVariables = [var.replace('_',' ') for var in meanVariables]

# format_row = '{:>30}'+'{:>12}'*(len(means[0]))
# formated_table = format_row.format('', *means[0])
# cell_text = []
# for var, val in zip(meanVariables, meanValues):
#    formated_table += '\n'
#    formated_table += format_row.format(var, *val[1])

table_meanValues=(r'\begin{tabular}{ c | c | c | c} '
+r' {} & {} & {} & {} \\\hline '.format(' ', *meanValues[0][0])
+r' {} & {} & {} & {} \\\hline '.format(meanVariables[0], *meanValues[0][1])
+r' {} & {} & {} & {} \\\hline '.format(meanVariables[1], *meanValues[1][1])
+r' {} & {} & {} & {} \\\hline'.format(meanVariables[2], *meanValues[2][1])
+r' {} & {} & {} & {} '.format(meanVariables[3], *meanValues[3][1])
+r'\end{tabular}' )
################# Calc Mean Values ############################

def getSubKeyList(substring):
   return [ key for key in FVS_dict.keys() if substring in key ]


FVS_data[:, FVS_dict['current_rpm']] = FVS_data[:, FVS_dict['current_rpm']] / ControlOptions['rpmmax']
GT_data[:, GT_dict['current_rpm']] = GT_data[:, GT_dict['current_rpm']] / ControlOptions['rpmmax']
print( "FVS maximum compressor_pressure = {}".format(np.max(FVS_data[FVS_t_index, FVS_dict['compressor_pressure']])) )
print( "GT maximum compressor_pressure = {}".format(np.max(GT_data[GT_t_index, GT_dict['compressor_pressure']])) )

# calculate the carnot efficiency following Kleins SEC_Code
FVS_eff_carnot = 1.0 - 1.0/FVS_data[FVS_t_index, FVS_dict['compressor_temperature']]
GT_eff_carnot = 1.0 - 1.0/GT_data[GT_t_index, GT_dict['compressor_temperature']]

f, axs = plt.subplots(nrows=4, ncols=2, figsize=(23/2,20/2) )
f.suptitle('Control Station' )

# add Latex table with mean values to plot
axs[0,0].text(1.5, 1.5, table_meanValues, ha="right", va="top", transform=axs[0,0].transAxes)


plotKeyList = ['mass_flow', 'power', 'pressure', 'temperature', 'rpm', 'fuel_consumption_rate', 'efficiency']#, 'SEC_Mode']
plotKeyList = ['mass_flow', 'power', 'pressure', 'temperature', 'rpm', 'fuel_consumption_rate', 'efficiency', 'efficiency']

data_units = ['', '', ' [bar]', ' [K]', ' [-]', '', '[-]', '[-]']
colors=['b', 'g', 'r', 'm']

for i, ax, subKey in zip(range(len(plotKeyList)), axs.flatten(), plotKeyList):
   subKeyList = getSubKeyList(subKey)
   for j, key in enumerate(subKeyList):
      if i>5: # don't destroy the labels from the last two plots
         lab = key
      else:
         lab = key.replace(subKey, '')
         lab = lab.replace('_', ' ')
      
      if (i==0):
         if ('comp' in lab):
            lab = lab.replace('compressor', 'comp. plenum')
         elif ('turb' in lab):
            lab = lab.replace('turbine', 'turb. plenum')
      
      if 'temperature' in subKey:
         FVS_data[FVS_t_index, FVS_dict[key]] *=T_ref
         GT_data[GT_t_index, GT_dict[key]] *=T_ref

      ax.plot( FVS_times[FVS_t_index], FVS_data[FVS_t_index, FVS_dict[key]], label=lab, color=colors[j] )
      ax.plot( GT_times[GT_t_index], GT_data[GT_t_index, GT_dict[key]], '--', zorder=10, color=colors[j] )

      if i>3:
         pass
      else:
         ax.legend()
      
      # if 'mass_flow' in subKey:
      #    ax.set(ylim=(0,1.2))

      if 'efficiency' in subKey:
         ax.plot(FVS_times[FVS_t_index], FVS_eff_carnot, label='Carnot efficiency', color=colors[1])
         ax.plot(GT_times[GT_t_index], GT_eff_carnot, '--', zorder=10, color=colors[1])
         if i==6:
            ax.legend()
            ax.set(ylim=(0,1.0))
      
      ylab = subKey.replace('_', ' ')
      if 'rpm' in ylab:
         ylab='rpm / rpmmax'
      ylab += data_units[i]
         
      ax.set(xlabel='time', ylabel=ylab, xlim=(FVS_times[0], None))
      ax.grid(True)

ax = axs.flatten()

# FVS_turb_density = FVS_data[:, FVS_dict['turbine_pressure']] / FVS_data[:, FVS_dict['turbine_temperature']] * rho_ref
# FVS_comp_density = FVS_data[:, FVS_dict['compressor_pressure']] / FVS_data[:, FVS_dict['compressor_temperature']] * rho_ref

# GT_turb_density = GT_data[:, GT_dict['turbine_pressure']] / GT_data[:, GT_dict['turbine_temperature']] * rho_ref
# GT_comp_density = GT_data[:, GT_dict['compressor_pressure']] / GT_data[:, GT_dict['compressor_temperature']] * rho_ref

# ax[-1].plot(FVS_times[FVS_t_index], FVS_comp_density[FVS_t_index], label='compressor', color=colors[0])
# ax[-1].plot( GT_times[GT_t_index], GT_comp_density[GT_t_index], '--', zorder=10, color=colors[0] )
# ax[-1].plot(FVS_times[FVS_t_index], FVS_turb_density[FVS_t_index], label='turbine', color=colors[1])
# ax[-1].plot( GT_times[GT_t_index], GT_turb_density[GT_t_index], '--', zorder=10, color=colors[1] )
# ax[-1].set(xlabel='time', ylabel='density '+r'$\left[ kg/m^3 \right]$', xlim=(FVS_times[0], None))
# ax[-1].legend()


f.subplots_adjust(hspace=0.3, wspace=0.3)
f.savefig('{}/control_state_FVS_GT.png'.format(output_path), bbox_inches='tight')

f.clear()
plt.close(f)
