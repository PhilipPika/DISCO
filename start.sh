#!/bin/zsh
emulate -LR zsh # reset zsh options
ulimit -S -n 2048
eval "$(/usr/local/anaconda3/bin/conda shell.zsh hook)";conda activate p37
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Mess: Exporting Dir for model"
export DGNM_ROOT=/Users/pippo/Desktop/DISCO/A_source_code
export DGNM_GENERALCODE=/Users/pippo/Desktop/DISCO/A_source_code/generalcode/trunk
export DGNM_USER=carbon

DIR=~/Desktop/DISCO
TOCODE=~/Desktop/DISCO/A_source_code/core

echo "Mess: Switching to Core Directory"
pwd
cd $TOCODE
pwd

echo "Mess: Removing current results"
rm -r $DIR/OUT

echo "Mess: Start SpinUp"
python dgnm_main.py --endtime=1951 --maskid=38 --lspinup=1 --inifile ../ini/cmd_m_100yrs_SpinUp_bio.ini
cp $DIR/OUT/bio/pkl/start1951.000.pkl $DIR/A_source_code/carbon/startups/start1951.000.pkl
echo "Note: First iteration done"

python dgnm_main.py --endtime=1960 --maskid=38 --lspinup=0 --inifile ../ini/cmd_m_100yrs_SpinUp_bio.ini
cp $DIR/OUT/bio/pkl/start1960.000.pkl $DIR/A_source_code/carbon/startups/start1951.000.pkl

echo "Mess: Start real run"
python dgnm_main.py --endtime=1970 --maskid=38 --lspinup=0 --inifile ../ini/cmd_m_50yrs_bio_def.ini
#
echo "Mess: Start Outout conversion"
cd $DIR/OUT
python ../A_source_code/carbon/code/output_conversion.py bio/pkl/ NETCDF
echo "Note: Outout conversion done"
