#!/bin/bash
echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"
echo "== Scratch dir. : ${TMPDIR}"

export PATH="/cm/shared/package/miniconda3/bin:$PATH"
source ~/.bashrc
ml miniconda/3

source activate /scistor/firelab/ppa760/.conda/envs/pipi37
conda activate /scistor/firelab/ppa760/.conda/envs/pipi37

## Define directories
DIR=${SLURM_SUBMIT_DIR}/DISCO
TOCODE=$DIR/A_source_code
export DGNM_ROOT=$TOCODE
export DGNM_GENERALCODE=$TOCODE/generalcode/trunk
export DGNM_USER=carbon

Addon=test_test_test
INI_file=species_bio_101 && echo "MESS: $INI_file"
StartUpEnd=1903 # Only StartUp
SimuEnd=1997 #Remember to also change the file name in params.ini
cmd_file=cmd_m_def.ini && echo "MESS: $cmd_file"

echo " "; echo "MESS: Removing current results" && rm -r $DIR/${Addon}_${INI_file}
echo "$DIR/${Addon}_${INI_file}";echo " "
########################################################################################
echo " ";cd $TOCODE/core && echo "MESS: Switching to Core Directory";echo " ";
python dgnm_main.py --lspinup=1 --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/
cp $DIR/${Addon}_${INI_file}/bio/pkl/start${StartUpEnd}.000.pkl $DIR/A_source_code/carbon/startups/start1901.000.pkl
python dgnm_main.py --lspinup=0 --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/
cp $DIR/${Addon}_${INI_file}/bio/pkl/start${StartUpEnd}.000.pkl $DIR/A_source_code/carbon/startups/start1901.000.pkl
#python dgnm_main.py --lspinup=0 --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/
#cp $DIR/${Addon}_${INI_file}/bio/pkl/start${SimuEnd}.000.pkl $DIR/A_source_code/carbon/startups/start1911.000.pkl

echo " ";cd $TOCODE/carbon/code/ && echo "MESS: START output conversion";echo " "
python output_conversion.py $DIR/${Addon}_${INI_file}/bio/pkl/ NETCDF

cd $TOCODE/core && echo "MESS: START aggragate TS"
python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=7 --mask_bool_operator=EQ
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=8 --mask_bool_operator=EQ
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=9 --mask_bool_operator=EQ
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=13 --mask_bool_operator=EQ
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=28 --mask_bool_operator=EQ
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=38 --mask_bool_operator=EQ
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=55 --mask_bool_operator=EQ
#python ../carbon/code/aggregate_timeseries.py --inifile ../ini/${cmd_file} --species_ini=${INI_file}.ini --endtime=$SimuEnd --outputdir=../../${Addon}_${INI_file}/bio/pkl/ --maskid=220 --mask_bool_operator=EQ

echo "MESS: FINISHED ${Addon}_$INI_file"
#########################################################################################
