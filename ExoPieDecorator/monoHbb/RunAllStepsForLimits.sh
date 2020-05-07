workspace=monoHbb_WS.root
datacardsDir=datacards_monoHbb_2017
pythonmacro=RunLimits.py
model=2hdma
category=merged

root -l -b -q PrepareWS.C
cp $workspace $datacardsDir
python $pythonmacro -c --model $model --$category --region "SR TOPE TOPMU WE WMU ZEE ZMUMU" 
#python $pythonmacro -c --model $model --$category --region "SR TOPMU" 
combine -M FitDiagnostics --saveShapes datacards_monoHbb_2017/datacard_monoHbb2017_B_Merged_sp_0p35_tb_1p0_mXd_10_mA_1000_ma_150.txt
python diffNuisances.py fitDiagnostics.root --abs --all


## long command for fit Diagnostics 
combine -M FitDiagnostics datacards_monoHbb_2017/datacard_monoHbb2017_B_Merged_sp_0p35_tb_1p0_mXd_10_mA_1000_ma_150.txt  --saveShapes --saveWithUncertainties --saveNormalizations --X-rtd MINIMIZER_analytic --rMin -100 -t -1 --expectSignal 0 -n _blind_wopos


combine -M AsymptoticLimits datacards_monoHbb_2017/datacard_monoHbb2017_B_Merged_sp_0p35_tb_1p0_mXd_10_mA_1000_ma_150_bkgsum.txt  --noFitAsimov  -v 3 --rMin 1e-07 --rMax 30  
