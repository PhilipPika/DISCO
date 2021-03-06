#!/bin/zsh
emulate -LR zsh # reset zsh options
ulimit -S -n 2048
eval "$(/usr/local/anaconda3/bin/conda shell.zsh hook)";conda activate p37
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "MESS: Exporting Dir for model"

DIR=/Users/pippo/Documents/SurfDrive/Research/Projects/DISCO_project/MyVersion/DISCO
TOCODE=$DIR/A_source_code

export DGNM_ROOT=$TOCODE
export DGNM_GENERALCODE=$TOCODE/generalcode/trunk
export DGNM_USER=carbon








Addon=test_2
INI_file=species_bio_101 && echo "MESS: $INI_file"
StartUpEnd=1999 # Only StartUp
SimuEnd=1999 #Remember to also change the file name in params.ini
cmd_file=cmd_m_def.ini && echo "MESS: $cmd_file"

#echo " "; echo "MESS: Removing current results" && rm -r $DIR/${Addon}_${INI_file}
echo "$DIR/${Addon}_${INI_file}";echo " "
#########################################################################################
echo " ";cd $TOCODE/core && echo "MESS: Switching to Core Directory";echo " ";
python dgnm_main.py --lspinup=1 --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$StartUpEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/
cp $DIR/${Addon}_${INI_file}/bio/pkl/start${StartUpEnd}.000.pkl $DIR/A_source_code/carbon/startups/start1901.000.pkl
python dgnm_main.py --lspinup=0 --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$StartUpEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/
cp $DIR/${Addon}_${INI_file}/bio/pkl/start${StartUpEnd}.000.pkl $DIR/A_source_code/carbon/startups/start1901.000.pkl
python dgnm_main.py --lspinup=0 --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/
cp $DIR/${Addon}_${INI_file}/bio/pkl/start${SimuEnd}.000.pkl $DIR/A_source_code/carbon/startups/start1901.000.pkl

echo " ";cd $TOCODE/carbon/code/ && echo "MESS: START output conversion";echo " "
python output_conversion.py $DIR/${Addon}_${INI_file}/bio/pkl/ NETCDF

cd $TOCODE/core && echo "MESS: START aggragate TS"
python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=7 --mask_bool_operator=EQ

echo "MESS: FINISHED ${Addon}_$INI_file"
#########################################################################################
