#!/bin/zsh
emulate -LR zsh # reset zsh options
ulimit -S -n 4096
eval "$(/usr/local/anaconda3/bin/conda shell.zsh hook)";conda activate p37
#export DGNM_USER=carbon;
#export DGNM_GENERALCODE=/work/philpika/DISCO_Test/A_source_code/generalcode/trunk;
#export DGNM_ROOT=/work/philpika/DISCO_Test/A_source_code
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

export DGNM_ROOT=/Users/pippo/Desktop/CARBON_DISC/A_source_code
export DGNM_GENERALCODE=/Users/pippo/Desktop/CARBON_DISC/A_source_code/generalcode/trunk
export DGNM_USER=carbon

DIR=~/Desktop/DISCO
TOCODE=~/Desktop/DISCO/A_source_code/core
rm -r $DIR/OUT
pwd
cd $TOCODE
pwd

#python dgnm_main.py --endtime=1951 --maskid=38 --lspinup=1 --inifile ../ini/cmd_m_50yrs_bio_def.ini
#cp $DIR/OUT/bio/pkl/start1951.000.pkl $DIR/A_source_code/carbon/startups/start1951.000.pkl
#echo "Note: First iteration done"

#for idx in {1..2};do
python dgnm_main.py --endtime=1951 --maskid=38 --lspinup=0 --inifile ../ini/cmd_m_50yrs_bio_def.ini
#cp $DIR/OUT/bio/pkl/start1951.000.pkl $DIR/A_source_code/carbon/startups/start1951.000.pkl
##python $TOCODE/../carbon/startups/ResetsYear.py
#echo "Note: Loop $idx done"
#done
#
#python dgnm_main.py --endtime=1951 --maskid=38 --lspinup=0 --inifile ../ini/cmd_m_50yrs_bio_def.ini
#
#cd $DIR/OUT
#python ../A_source_code/carbon/code/output_conversion.py bio/pkl/ NETCDF
#echo "Note: Outout conversion done"
