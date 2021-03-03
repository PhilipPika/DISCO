#!/bin/zsh
emulate -LR zsh # reset zsh options
ulimit -S -n 2048
eval "$(/usr/local/anaconda3/bin/conda shell.zsh hook)";conda activate p37
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Mess: Exporting Dir for model"

DIR=/Users/pippo/Documents/SurfDrive/Research/Projects/DISCO_project/MyVersion/DISCO
TOCODE=$DIR/A_source_code

export DGNM_ROOT=$TOCODE
export DGNM_GENERALCODE=$TOCODE/generalcode/trunk
export DGNM_USER=carbon

echo "Mess: Switching to Core Directory"
cd $TOCODE/core

echo "Mess: Removing current results"
rm -r $DIR/OUT
Spinstart=1988
SpinUpEnd=1990
SimuEnd=1999 #Remember to also change the file name in params.ini

echo "Mess: Start SpinUp"
nice -n 19 python dgnm_main.py --lspinup=1 --inifile ../ini/Ccmd_m_50yrs_bio_def.ini --endtime=$SpinUpEnd
cp $DIR/OUT/bio/pkl/start$SpinUpEnd.000.pkl $DIR/A_source_code/carbon/startups/start$Spinstart.000.pkl

cp $DIR/OUT $DIR/OUT_nonSS

echo "Note: First iteration done"
nice -n 19 python dgnm_main.py  --lspinup=0 --inifile ../ini/Ccmd_m_50yrs_bio_def.ini

echo "Mess: Start Outout conversion"
cd $DIR/OUT
python ../A_source_code/carbon/code/output_conversion.py bio/pkl/ NETCDF
echo "Note: Outout conversion done"

echo "Mess: Start aggregate time series"
cd $TOCODE/core
python ../carbon/code/aggregate_timeseries.py --inifile ../ini/Ccmd_m_50yrs_bio_def.ini --endtime=$SimuEnd
echo "Mess: Aggregate time series done"
