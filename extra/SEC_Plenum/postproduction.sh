CaseName=$1 # get path from cli
NCPUs=8 # for parallel plotting

time python3 PlotTube_HDF5.py $CaseName
time python3 PlotControlState_HDF5.py $CaseName
# time python3 PlotPlenum_HDF5.py $CaseName --parallel=$NCPUs # dauert lange
time python3 PlotPlenumStats.py $CaseName