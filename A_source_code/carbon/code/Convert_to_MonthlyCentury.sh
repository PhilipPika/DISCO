#!/bin/zsh
pp37 

cd /Users/pippo/Documents/SurfDrive/Research/Projects/DISCO_project/MyVersion/DISCO/A_source_code/carbon/code
python convert_nc.py 45_to_101 ../../../../Monthly_B_model_input/c_delivery/soil_loss_OC_mon.nc
python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/soil_loss_OC_mon_101.nc

python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/C_gnpp_headwaters_101.nc
python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/ALK_TOTAL_flux_yr_101.nc
python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/doc_sewage_101.nc
python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/soil_loss_TSS_yr_101.nc
python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/phyto_ini_101.nc
python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/C_NPP_HarvestCORR_floodplain_101.nc
python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/DOC_outflux_for_total_gridcell_area_101.nc
python convert_nc.py 101_to_1200 ../../../../Monthly_B_model_input/c_delivery/doc_sewage_101.nc
cd -
