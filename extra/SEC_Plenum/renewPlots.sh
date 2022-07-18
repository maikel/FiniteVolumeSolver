
file=Plenum_Pressure_stats.dat
for dir in ../../build_2D-Release/sec_parameterstudy/*; 
do
    if [[ -d "$dir" && ! -L "$dir" ]]; then
        echo "$dir"
        # if [ ! -f "$dir/Visualization/$file" ]; then
        #     echo "$dir"
        #     # python3 ComputePlenumStats.py "$dir" "SEC_Plenum_Arrhenius_oneTube.py"
        #     #ls "$dir/Visualization/" -l
        # fi
        
        python3 ComputePlenumStats.py "$dir" SEC_Plenum_Arrhenius_oneTube.py
        python3 PlotPlenumStats.py "$dir" SEC_Plenum_Arrhenius_oneTube.py
	python3 PlotTube_HDF5.py "$dir" SEC_Plenum_Arrhenius_oneTube.py
        python3 PlotControlState_HDF5.py "$dir" SEC_Plenum_Arrhenius_oneTube.py
    fi 
done
