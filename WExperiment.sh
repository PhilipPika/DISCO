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

cd $TOCODE/core; echo "MESS: Switching to Core Directory"

echo "MESS: Removing current results"
rm -r $DIR/OUT
Spinstart=1970
SpinUpEnd=1972
SimuEnd=1972 #Remember to also change the file name in params.ini
########################################################################################
INI_file=mon_bio_1200.ini
#INI_file=mon_bio_1200_2tPOCspecies.ini
#INI_file=mon_bio_1200_2tPOCspecies_transPOC2DOC.ini
echo "MESS: $INI_file"
#echo "nice -n 19 python dgnm_main.py --lspinup=1 --inifile ../ini/Ccmd_m_50yrs_bio_def.ini --endtime=$SpinUpEnd --species_ini=$INI_file"
nice -n 19 python dgnm_main.py --lspinup=1 --inifile ../ini/Ccmd_m_50yrs_bio_def.ini --endtime=$SpinUpEnd --species_ini=$INI_file
cp $DIR/OUT/bio/pkl/start$SpinUpEnd.000.pkl $DIR/A_source_code/carbon/startups/start$Spinstart.000.pkl
cp $DIR/OUT/bio/pkl/start$SpinUpEnd.000.pkl $DIR/A_source_code/carbon/startups/start$Spinstart.000.$INI_file.pkl
#echo "MESS: START Actual run"
nice -n 19 python dgnm_main.py  --lspinup=0 --inifile ../ini/Ccmd_m_50yrs_bio_def.ini --species_ini=$INI_file
echo "MESS: START output conversion"
cd $TOCODE/carbon/code/
python output_conversion.py $DIR/OUT/bio/pkl/ NETCDF
echo "MESS: START aggragate TS"
cd $TOCODE/core
python ../carbon/code/aggregate_timeseries.py --inifile ../ini/Ccmd_m_50yrs_bio_def.ini --endtime=$SimuEnd --species_ini=$INI_file
#mv $DIR/OUT $DIR/$INI_file
echo "MESS: FINISHED $INI_file"
#########################################################################################
