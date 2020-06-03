# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/dgnm_main.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Set the general path for the own python modules
import general_path

# Import general python modules
import optparse
import os
import pickle
import sys
import time
import traceback

# Import own general modules
import ascraster
from error import *
from iround import *
import my_sys

# Import local modules.
import define_subgrid_streamorder
import dgnm_river
import file_locking
import general_func
import general_startup
import make_cluster_cmd
import make_index_species
import make_photo
import read_parameter
import create_outfiles
import make_slurm_script

def run_dgnm(args):
    '''
    Run DGNM model main routine
    args: list of arguments passed trough command-line or other script
    '''

    # Parse command-line arguments and set parameters for script
    # Startup logging and runtime start
    params,log,s = general_startup.general_startup(args)

    # Read mask of the input grids.
    # The mask takes care to do the calculations on an efficient way.
    # Only the grid cells which are in the mask are calculated and stored.
    if (params.lmask):

        # Global
        #mask = ascraster.create_mask(params.file_mask, 0,'GT',numtype=float)
        # General formulation in the ini file.
        mask = ascraster.create_mask(params.file_mask, params.maskid, params.mask_bool_operator,numtype=int)
      
        if (params.ldebug): 
            print("MASK: ",params.file_mask  + " has been read.")
        log.write_and_print(s.interval("Reading mask."))
    else:    
        mask = None
        log.write_and_print(s.interval("No mask is used for this simulation."))
    
    # Read species file, conservative variables, sources and processes file
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)

    # Add all params_local to params.
    for iname in params_local.names:
        setattr(params,iname,params_local.get_val(iname))
    del params_local

    # Make an indexnumber in the species list with the name of the species.
    # For example: params.io2   = general_func.find_spec(species,"O2")   
    make_index_species.make_index_species(params,species,proc)

    # Make the photosythesis for each species which needs it.
    # Calculate this for each day and for just 1 year.
    # Calculation is only done if no photosynthesis input file is provided.
    make_photo.make_photosyn(params,species)  
    
    # Define characteristics of low Strahler orders if constant.
    if (not params.lstrahlergrids): #LV 28-02-2018
        define_subgrid_streamorder.define_subgrid_streamorder_constant(params)

    # Delete all the locking files on the output directory
    file_locking.clean_locks(params)
    
    # Read river basin map
    basin = ascraster.Asciigrid(ascii_file=params.basin,mask=mask,numtype=int)
    if (params.ldebug): 
        print(params.basin  + " has been read.")

    # Create the information of the mask to an output file to be able to read the output files.
    fp_out = open(os.path.join(params.outputdir,"mask.pkl"),"wb")
    pickle.dump(basin,fp_out, -1)
    pickle.dump(params,fp_out, -1)
    fp_out.close()

    # Make a list of all the riverbasins that have to be calculated.
    list_riverid = list(set(basin.values))
    list_riverid.sort()

    # Creation of pickle file with fixed initial concentrations (read in ini file).
    if  (params.lspinup == -1):
        general_func.make_header(os.path.join(params.pkl_outputdir,"start_fixed_conc.pkl"),species,params.starttime)

    elif  (params.lspinup == 0):
        # Check the startup file
        general_func.check_header(params.startup_file,species,params.starttime)
 
    elif  (params.lspinup == 1):
        # Make empty output files where all the processes can add their information.
        create_outfiles.create_outfiles(params,species,proc,params.starttime)

    # Make all output files.
    year = params.starttime + params.outputtime
    timeperiod = -1
    while (year < params.endtime):
        timeperiod += 1
        # Make all restart files.
        lrestart = ((timeperiod+1)%iround(params.restartnumber) == 0)
        create_outfiles.create_outfiles(params,species,proc,year,lrestart=lrestart)
        year += params.outputtime

    # Make output files for last time
    year = params.endtime
    create_outfiles.create_outfiles(params,species,proc,year)

    # Determine the number of cells per riverbasin, to avoid too small portions.
    riverid_num = {}
    for icell in range(basin.length):
        id = basin.get_data(icell)
        if (id != None):
            try:
                riverid_num[id] += 1
            except KeyError:
                riverid_num[id] = 1

    # Store the objects that are needed in the dgnm_river.
    # Write these objects to file, so dgnm_river can read this file.
    # Here the general data is accessible for all subprocesses.
    fp_out = open(os.path.join(params.outputdir,"general_objects.pkl"),"wb")
    pickle.dump(params,fp_out, -1)
    pickle.dump(species,fp_out, -1)
    pickle.dump(proc,fp_out, -1)
    pickle.dump(sources,fp_out, -1)
    fp_out.close()


    # Prepare for multitasking.
    if (params.llocal):

        # Use the local processing power and use the multiprocessing module of python.
        import multiprocessing as mp

        # Define a lock to regulate the access of shared objects over all the processes.
        lock = mp.Lock()

        # Set the number of cpu that will be used.
        number_of_cpu = min(mp.cpu_count(),params.ncpu)

        # Create list with jobs.       
        jobs= []

    # Start with submitting the jobs
    item = 0
    jobnum = 1
    while (item < len(list_riverid)):
        try:
            total_cells = riverid_num[list_riverid[item]]
        except KeyError:
            item += 1
            continue
        numbers = [list_riverid[item]]

        while (total_cells < params.minimal_number_of_cells):
            item += 1
            if (item > len(list_riverid)-1):
                # Last riverbasin, so no more cells.
                break
            else:
                total_cells += riverid_num[list_riverid[item]]
                numbers.append(list_riverid[item])
        if (params.llocal):
            # Define the process with riverid "num", return value p is process "handler"
            p=mp.Process(target= dgnm_river.calculate, args=(numbers,lock,params,species,proc,sources))
            # Add the process handler to the jobs list.
            jobs.append(p)
            # Start the process
            general_func.start_process(p,jobs,number_of_cpu)

        else:
            # Distribute the work over the nodes of the EEJIT cluster with sbatch
            # Make command string to execute on the cluster
            cmd = make_cluster_cmd.make_cluster_cmd(params,jobnum,numbers,total_cells)
            # Submit the job.
            my_sys.my_system(cmd)
            time.sleep(2)
            
        item += 1
        jobnum += 1 

    # Wait until all jobs are finished.
    if (params.llocal):    
        for p in jobs:
            p.join()
        print("All processes are ready.")
    
    log.write_and_print(s.total("Total run"))    
    del log

if __name__ == "__main__":
  
    try:
        starttime_main = time.time()

        run_dgnm(sys.argv)

        # end timer  
        endtime_main = time.time()
        print('Total simulation time (in s):  ' + str(endtime_main-starttime_main))
    except MyError as val:
        val.write()

    except (optparse.OptionValueError,
            ascraster.ASCIIGridError,
            Exception) as val:
        print(str(val))
        print(str(sys.exc_info()[0]))
        print(traceback.print_exc())
    except:
        print("***** ERROR ******")
        print("dgnm_main.py failed.")
        print(str(sys.exc_info()[0]))
        print(traceback.print_exc())
