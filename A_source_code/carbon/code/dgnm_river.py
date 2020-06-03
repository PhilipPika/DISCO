# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Set the general path for the own python modules
try:
    # Try to load the general_path module to set the python path for own python modules.
    import general_path
except ImportError:
    import os
    import sys
    root = os.getenv("DGNM_ROOT")
    if (root == None):
        print("***** ERROR ******")
        print("Environment parameter DGNM_ROOT is not set.")
        sys.exit(1)
    if (not os.path.isdir(root)):
        print("***** ERROR ******")
        print("Environment parameter DGNM_ROOT is not set correctly.")
        print("Environment parameter DGNM_ROOT found: ",root)
        sys.exit(1)

    # Put core directory in the sys.path.
    path = os.path.join(root,"core")
    if (os.path.exists(path)): 
        sys.path.insert(0, path)

    # Try it again to load the general_path module.
    import general_path


# Import Python modules
import optparse
import os
import pickle
import random
import sys
import time
import traceback

# Import own general modules
import accuflux
import ascraster
from error import *
from iround import *
import my_sys
import pointer_class

# Import local modules.
import dgnm_cell
import file_locking
import general_func
import get_dynamic_gridinfo
import interpolate_list
import make_outlakes_dict
import define_subgrid_streamorder
import vol_depth_vel
import directory

def calculate(iriv,lock,params,species,proc,sources):
    '''
    Main function to calculate a riverbasin or a set of riverbasins.
    '''
    # Get the river id's which are being calculated.
    if (type(iriv) == list):
        rivlist = iriv
    else:
        rivlist = [iriv]
    print("Rivernumbers: ", rivlist)

    # Make a dictionary with river ids for each cell in the mask.
    rivernum = {}
    # Make mask for this collection of basins.
    rivmask = []
    cellnumber = 0

    # Read grid file with all basin ids without mask.
    basin = ascraster.Asciigrid(ascii_file=params.basin,numtype=int)

    # Read grid file with all basin ids.
    for icell in range(basin.length):
        id = basin.get_data(icell)
        if (id != None):
            if (iround(id) in rivlist):
                rivmask.append(icell)
                rivernum[cellnumber] = iround(id)
                cellnumber += 1
        
    # Set jobid for this collection of riverbasins.
    jobid = rivlist[0]
    print("Reading mask for riverids with jobid",jobid)

    # Read ldd map 
    lddmap = ascraster.Asciigrid(ascii_file=params.ldd,mask = rivmask,numtype=int)
    print("Reading ldd for riverid ",jobid)

    # Make a pointer for the calculation order
    pointer1 = pointer_class.make_pointer(params, rivmask=rivmask)

    # Set temp output directory for this river basin.
    # First check the temp directory
    #params.tmp_dir = directory.ensure(os.path.join(params.outputdir, params.tmp_dir))

    #unique_jobid = str(int(random.random()*1000))+"_"+str(jobid)
    tmpdir = os.path.join(params.tmp_dir, str(jobid))

    # Create temp output directory for this riverbasin
    #if os.path.exists(tmpdir):
        # Remove old directory.
    #    my_sys.my_rmdir(tmpdir)
    # Make new temporary output directory
    os.mkdir(tmpdir)

    print(tmpdir)

    # Create photosynthese files for each species which needs photosynthese.
    if (len(params.totalphotindex) > 0):  
        photo_list = []
        for ispec in params.totalphotindex:
            specname = species[ispec].get_name()
            fp = open(os.path.join(params.outputdir,\
                            os.path.splitext(os.path.basename(params.photosyn_file))[0]\
                             + "_" + specname + ".pkl"),"rb")

            photo_list.append(list(general_func.pickleLoader(fp)))
            fp.close()
            
        for item in range(len(pointer1)):
            icell = pointer1[item].get_index()
            # Get the latitude from the cell
            lon,lat = basin.get_coordin_from_index(icell)
            ind = iround((90.0 - lat+0.5*basin.cellsize)/basin.cellsize)
            fp_out = open(os.path.join(tmpdir,"photo_"+str(pointer1[item].get_local_index())+".pkl"),"wb")
            for i in range(len(params.totalphotindex)):
                pickle.dump(photo_list[i][ind][1:],fp_out,-1)
            fp_out.close()

    # Read the constant channel_width of PCRGLOBWB. # added in volume_depth.... calculations
    channel_width_grid = ascraster.Asciigrid(ascii_file=params.channel_width,\
                                             mask=rivmask,numtype=float)
            
    # Read the slopes of the cells
    slopemap = ascraster.Asciigrid(ascii_file=params.slope,\
                                   mask=rivmask,numtype=float)

    if params.lfloodplains==1:
        norders = params.norder+1
    else:
        norders = params.norder

    #outlakes_dict = make_outlakes_dict.calculate(params,pointer1)
            
    # Get all data from the input files with multiprocessing.
    if (params.ncpu_read_files > 1):

        # Make preparations for multiprocessing this part of the code.
        # Use the local processing power and use the multiprocessing module of python.
        import multiprocessing as mp

        # Set the number of cpu that will be used.
        number_of_cpu = min(mp.cpu_count(),params.ncpu_read_files)

        # Make a list of jobs which are active               
        jobs= []
            
        outlakes_dict = make_outlakes_dict.calculate(params,pointer1)

        # Calculate volume and depth.
        #p=mp.Process(target= volume_depth_netcdf.calculate, args=(params,rivmask,tmpdir,pointer1,str(jobid),True))
        #p=mp.Process(target= volumes_depths.calculate, args=(params,rivmask,tmpdir,pointer1,outlakes_dict,str(jobid),True))
        p=mp.Process(target= vol_depth_vel.calculate, args=(params,rivmask,tmpdir,pointer1,outlakes_dict,str(jobid),True)) #LV 14-02-2018
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)

        # Read lake/reservoir id's for the whole period and put this time dependent object.
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo,\
                     args=(params,pointer1,os.path.join(params.water_inputdir,params.lakeid),params.lakeid_varname,None,\
                           os.path.join(tmpdir,"lakeid_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)
        
        # Read lake/reservoir outflow id's for the whole period and put this time dependent object.
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo,\
                     args=(params,pointer1,os.path.join(params.water_inputdir,params.endo_lakes),params.endo_lakes_varname,None,\
                           os.path.join(tmpdir,"endo_lakes_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)
            
        # Read endoreic lake/reservoir id's for the whole period and put this time dependent object.
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo,\
                     args=(params,pointer1,os.path.join(params.water_inputdir,params.outlakeid),params.outlakeid_varname,None,\
                           os.path.join(tmpdir,"outlakeid_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)
            
        # Read the discharge for the whole period and put this time dependent object km3/yr.
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo,\
                     args=(params,pointer1,os.path.join(params.water_inputdir,params.discharge),params.discharge_varname,None,\
                           os.path.join(tmpdir,"discharge_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)

        # Read the surface water area of lakes and resrvoirs m2
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo,\
                     args=(params,pointer1,os.path.join(params.water_inputdir,params.water_area),params.water_area_varname,None,\
                           os.path.join(tmpdir,"water_area_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)
            
        # Read the runoff for the whole period and put this time dependent object mm/yr.
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo,\
                     args=(params,pointer1,os.path.join(params.water_inputdir,params.pnet),params.pnet_varname,None,\
                           os.path.join(tmpdir,"runoff_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)

        # Read the temperature for the whole period and put this time dependent object in degree Celcius.
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo,\
                     args=(params,pointer1,os.path.join(params.water_inputdir,params.temperature),params.temperature_varname,None,\
                           os.path.join(tmpdir,"temperature_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)

        # Read the global radiation for the whole period and put this time dependent object in W/m2.
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo, args=(params,pointer1,\
                             os.path.join(params.water_inputdir,params.global_radiation),params.global_radiation_varname,None,\
                             os.path.join(tmpdir,"global_radiation_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)

        # Read the low vegetation fraction for the whole period and put this time dependent object in dimensionless
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo, args=(params,pointer1,\
                             os.path.join(params.water_inputdir, '..', 'other', params.low_veg_fr),params.low_veg_fr_varname,None,\
                             os.path.join(tmpdir,"low_veg_fr_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)

        # Read the high vegetation fraction for the whole period and put this time dependent object in dimensionless
        p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo, args=(params,pointer1,\
                             os.path.join(params.water_inputdir, '..', 'other', params.high_veg_fr),params.high_veg_fr_varname,None,\
                             os.path.join(tmpdir,"high_veg_fr_"+str(jobid)+".pkl")))
        # Add the process handler to the jobs list.
        jobs.append(p)
        # Start the process
        general_func.start_process(p,jobs,number_of_cpu)

        # Read loading for all the sources
        for isrc in range(len(sources)):
            filename = os.path.join(params.load_inputdir,sources[isrc].get_val('name')+'.nc')
            temp_distrib = None
            if ('temporal' in sources[isrc].get_attrib()):
                temp_distrib = sources[isrc].get_val('temporal')
            p=mp.Process(target= get_dynamic_gridinfo.get_dynamic_gridinfo,\
                         args=(params,pointer1,filename,sources[isrc].get_val('name'),temp_distrib,\
                               os.path.join(tmpdir,sources[isrc].get_val('name')+"_"+str(jobid)+".pkl")))
            # Add the process handler to the jobs list.
            jobs.append(p)
            # Start the process
            general_func.start_process(p,jobs,number_of_cpu)

        # Wait until all jobs are finished.
        for p in jobs:
            p.join()

        # Get all the content of the files which dumped to file.
        # LV 06-02-2017: added volume and depth of water in floodplains
        #volume_cells = general_func.get_content_rm(os.path.join(tmpdir,"volume_"+str(jobid)+".pkl"))
        #depth_cells = general_func.get_content_rm(os.path.join(tmpdir,"depth_"+str(jobid)+".pkl"))
        volume_cells = general_func.get_content_rm(os.path.join(tmpdir,"volume_c_"+str(jobid)+".pkl"))
        depth_cells = general_func.get_content_rm(os.path.join(tmpdir,"depth_c_"+str(jobid)+".pkl"))
        volume_fp_cells = general_func.get_content_rm(os.path.join(tmpdir,"volume_fp_"+str(jobid)+".pkl"))
        depth_fp_cells = general_func.get_content_rm(os.path.join(tmpdir,"depth_fp_"+str(jobid)+".pkl"))
        vel_cells = general_func.get_content_rm(os.path.join(tmpdir,"velocity_"+str(jobid)+".pkl")) #LV 14-02-2018
        water_area_cells = general_func.get_content_rm(os.path.join(tmpdir,"water_area_"+str(jobid)+".pkl")) #LV 29-11-2017
        runoff_cells = general_func.get_content_rm(os.path.join(tmpdir,"runoff_"+str(jobid)+".pkl"))
        lakeid_cells = general_func.get_content_rm(os.path.join(tmpdir,"lakeid_"+str(jobid)+".pkl"))
        outlakeid_cells = general_func.get_content_rm(os.path.join(tmpdir,"outlakeid_"+str(jobid)+".pkl"))
        endo_lakes_cells = general_func.get_content_rm(os.path.join(tmpdir,"endo_lakes_"+str(jobid)+".pkl"))
        discharge_cells = general_func.get_content_rm(os.path.join(tmpdir,"discharge_"+str(jobid)+".pkl"))
        temperature_cells = general_func.get_content_rm(os.path.join(tmpdir,"temperature_"+str(jobid)+".pkl"))
        globrad_cells = general_func.get_content(os.path.join(tmpdir,"globrad_"+str(jobid)+".pkl"))
        low_veg_fr_cells = general_func.get_content(os.path.join(tmpdir,"low_veg_fr_"+str(jobid)+".pkl"))
        high_veg_fr_cells = general_func.get_content(os.path.join(tmpdir,"high_veg_fr_"+str(jobid)+".pkl"))
            
        srcload = []
        for isrc in range(len(sources)):
            filename = os.path.join(tmpdir,sources[isrc].get_val('name')+"_"+str(jobid)+".pkl")
            srcload.append(general_func.get_content(filename))                
                
    else:
            
        outlakes_dict = make_outlakes_dict.calculate(params,pointer1)

        # Calculate volume and depth.
        #volume_cells,depth_cells = volume_depth_netcdf.calculate(params,rivmask,tmpdir,pointer1,str(jobid),False)
        #volume_cells,depth_cells, volume_fp_cells, depth_fp_cells = volumes_depths.calculate(params,rivmask,tmpdir,pointer1,outlakes_dict,str(jobid),False)
        volume_cells,depth_cells, volume_fp_cells, depth_fp_cells, vel_cells = vol_depth_vel.calculate(params,rivmask,tmpdir,pointer1,outlakes_dict) #LV 14-02-2018
        #print "vel_cells[3]  " + str(vel_cells[3])
            
        # Read surface water area of lakes and reservoirs m2
        water_area_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir,params.water_area),params.water_area_varname) #LV 29-11-2017
            
        # Read the runoff for the whole period and put this time dependent object mm/yr.
        runoff_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir,params.pnet),params.pnet_varname)
            
        # Read the lake information for the whole period and put this time dependent object.
        lakeid_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir,params.lakeid),params.lakeid_varname)
        outlakeid_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir,params.outlakeid),params.outlakeid_varname)
        endo_lakes_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir,params.endo_lakes),params.endo_lakes_varname)

        # Read the discharge for the whole period and put this time dependent object km3/yr.
        discharge_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir,params.discharge),params.discharge_varname)

        # Read the temperature for the whole period and put this time dependent object in degree Celcius.
        temperature_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir,params.temperature),params.temperature_varname)

        # Read the global radiation for the whole period and put this time dependent object in degree W/m2
        globrad_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir,params.global_radiation),params.global_radiation_varname)

        low_veg_fr_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir, params.low_veg_fr),params.low_veg_fr_varname)

        high_veg_fr_cells = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(params.water_inputdir, params.high_veg_fr),params.high_veg_fr_varname)

        # Read loads from all sources and apply temporal pattern
        srcload = []
        for isrc in range(len(sources)):
            filename = os.path.join(params.load_inputdir,sources[isrc].get_val('name')+'.nc')
            temp_distrib = None
            if ('temporal' in sources[isrc].get_attrib()):
                temp_distrib = sources[isrc].get_val('temporal')
            srcload.append(get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,filename,\
                                                                     sources[isrc].get_val('name'),\
                                                                     temp_distrib=temp_distrib))

    # Create strahler argument list if gridded chracteristics
    args_strahler_cell = None
    if (params.lstrahlergrids):
        args_strahler = define_subgrid_streamorder.define_subgrid_streamorder_riverbasin(params,rivmask)

    # Make preparations for multiprocessing of riverbasin dynamic calculations.
    number_of_cpu = general_func.get_number_of_cores(params,len(pointer1))

    if (params.ncpu_riverbasin > 1 and number_of_cpu > 1):
        # Make preparations for multiprocessing this part of the code.
        # Use the local processing power and use the multiprocessing module of python.
        import multiprocessing as mp

        # Set the number of cpu that will be used.
        number_of_cpu = min(mp.cpu_count(),number_of_cpu+1)

        # Make a list of jobs which are active               
        jobs= []
     
        # Prev_order is needed to find out which cells can be run at the same time. Cells of the same order can run at the same time.
        prev_order = pointer1[0].get_iorder()
        
        # Make lock for write and reading loads.
        lock_cell = {}
    else:
        number_of_cpu = 1

    # Read start file and extract information for all cells in this riverbasin and
    # store all information of the cells in separate input files.
    # In case of lspinup == 1 there is no startup file, so we have to create one ourselves.
    if ((params.lspinup == -1) or (params.lspinup == 0)):
        lfound = False
        # Open output file for first cell of this riverbasin. 
        icell_prev = pointer1[0].get_index()
        local_cell_number = rivmask.index(icell_prev) 
        fp_out = open(os.path.join(tmpdir,str(local_cell_number)+"_-1.pkl"),"wb")

        # Make a dictionary which indicates whether a river is found or not.
        lfound = {}
        for i in rivlist:
            lfound[i] = False

        # Read the startup file of all the riverbasins.  
        with open(params.startup_file, 'rb') as fp:
            # Read header and skip this
            header_skip = pickle.load(fp)

            # Read header and skip this
            general_func.pickleLoader(fp)

            for line in general_func.pickleLoader(fp):
                # Find first line with the riverid iriv, then lfound will be true.
                if (line[0] in rivlist):
                    # Now the start of the riverbasin is found.
                    lfound[iround(line[0])] = True
                    icell = line[1]
                    if (icell != icell_prev):
                        # Now there is another cell found. So close the previous file and open a new one for the new cell.
                        fp_out.close()
                        local_cell_number = rivmask.index(icell) 
                        fp_out = open(os.path.join(tmpdir,str(local_cell_number)+"_-1.pkl"),"wb")

                    # Put starting year in the file by overwriting the cellnumber.
                    line[1] = params.starttime
                    # Write cell info (starttime, and info about all the specs to output file.
                    pickle.dump(line[1:],fp_out,-1)
                    icell_prev = icell

            fp_out.close()
 
        lerror = False
        for key in rivlist:
            if (not lfound[key]):
                print("Startup file has no information on riverbasin "+str(key))
                lerror = True
        if (lerror): 
            raise MyError("Startup file has no information for riverbasins.") 

    else:
        # We have to make an inital condition
        lsteady = True
        timeperiod = -1
        yearstart = params.starttime
        yearend   = min(yearstart + params.outputtime,params.endtime)
        #lakes_dict = {}
        #outlakes_dict = {}
        # Loop over all cells:    
        for item in range(len(pointer1)):
            #icell = pointer1[item][1]
            icell = pointer1[item].get_local_index()
            direction = lddmap.get_data(icell)
            next_cell = accuflux.goto_nextcell(icell,direction,lddmap.ncols,mask=rivmask)

            # Load is the time array of  all the sources       
            load = []
            for isrc in range(len(sources)):
                load.append(srcload[isrc][icell])      

            # Check whether cell is a lake/reservoir.
            lakeid = interpolate_list.stepwise(yearstart,lakeid_cells[icell],extrapol=1)[-1]
            if (lakeid > 0):
                outlakeid = interpolate_list.stepwise(yearstart,outlakeid_cells[icell],extrapol=1)[-1]
                endo_lake = interpolate_list.stepwise(yearstart,endo_lakes_cells[icell],extrapol=1)[-1]
                llake = True
                llakeout = (outlakeid  > 0.)
                lendolake = (endo_lake > 0.)
                #lakes_dict[icell] = iround(lakeid)
                #if (llakeout or lendolake):
                #    outlakes_dict[iround(lakeid)] = icell
            else:
                llake = False
                llakeout = False
                lendolake = False
        
            if (number_of_cpu > 1):

                try:
                    qq = lock_cell[next_cell]
                except KeyError:
                    lock_cell[next_cell] = mp.Lock()
            
                if (prev_order != pointer1[item].get_iorder()):
                    # We have to wait untill all processes are ready, before we can continue
                    for p in jobs:
                        p.join()
                    prev_order = pointer1[item].get_iorder()
                    jobs = []

                # Define the process with cellid "icell", return value p is process "handler"
                width = channel_width_grid.get_data(icell)
                slope = slopemap.get_data(icell)

                if (params.lstrahlergrids):
                    args_strahler_cell = args_strahler[icell]
                p=mp.Process(target= dgnm_cell.calculate_cell,\
                             args=(lsteady,lock_cell[next_cell],icell,params,species,proc,sources,\
                                   tmpdir,next_cell,yearstart,yearend,timeperiod,\
                                   runoff_cells[icell],load,temperature_cells[icell],\
                                   globrad_cells[icell], discharge_cells[icell],\
                                   volume_cells[icell],water_area_cells[icell],depth_cells[icell],width,slope,\
                                   volume_fp_cells[icell],depth_fp_cells[icell],vel_cells[icell],\
                                   low_veg_fr_cells[icell],high_veg_fr_cells[icell],llake,llakeout,lendolake,args_strahler_cell))


                # Add the process handler to the jobs list.
                jobs.append(p)
                # Start the process
                general_func.start_process(p,jobs,number_of_cpu)
                
            else:
                # Do the job one at the time.
                width = channel_width_grid.get_data(icell)
                slope = slopemap.get_data(icell)
 
                if (params.lstrahlergrids):
                    args_strahler_cell = args_strahler[icell]
                dgnm_cell.calculate_cell(lsteady,None,icell,params,species,proc,sources,\
                                              tmpdir,\
                                              next_cell,yearstart,yearend,timeperiod,\
                                              runoff_cells[icell],load,temperature_cells[icell],\
                                              globrad_cells[icell], discharge_cells[icell],\
                                              volume_cells[icell],water_area_cells[icell],depth_cells[icell],width,slope,\
                                              volume_fp_cells[icell],depth_fp_cells[icell],vel_cells[icell],\
                                              low_veg_fr_cells[icell],high_veg_fr_cells[icell],llake,llakeout,lendolake,args_strahler_cell=args_strahler_cell)
            
        if (number_of_cpu > 1):            
            # Wait until all jobs are finished.
            for p in jobs:
                p.join()

        # Put cell information of the last time step into the world map.
        data_block = []
        data_block1 = []
        outputlist = []
        startoutputlist = []
        conclist = []
        budlist = []
        argslist = []
        data_block2 = []
        data_block3 = []

        # Read the pickled data for individual cells and store in one list (per river basin)
        for item in range(len(pointer1)):
            icell = pointer1[item].get_local_index()

            # For amounts
            with open(os.path.join(tmpdir,str(icell)+"_"+str(timeperiod)+".pkl"),"rb") as fp:
                for line in general_func.pickleLoader(fp):
                    data_block.append(line)
            fp.close()
            # Store only the last output time of the cell.
            for i in range(len(data_block)-norders,len(data_block)):
                # Put riverid and cellnumber in front of the line and remove the time.
                linestart = [rivernum[icell],rivmask[icell]]
                linestart.extend(data_block[i][1:])
                startoutputlist.append(linestart)
                # To be able to print species amounts/concentrations for all orders
                if (params.allorders == 1):
                    outputlist.append(linestart)
                else:
                    if (i == (len(data_block)-1)):
                        outputlist.append(linestart)
            
            # For concentrations
            with open(os.path.join(tmpdir,"conc_" + str(icell)+"_"+str(timeperiod)+".pkl"),"rb") as fp:
                for line in general_func.pickleLoader(fp):
                    data_block1.append(line)
            fp.close()
            for i in range(len(data_block1)-norders,len(data_block1)):
                # Put riverid and cellnumber in front of the line and remove the time.
                linestart = [rivernum[icell],rivmask[icell]]
                linestart.extend(data_block1[i][1:])
                # To be able to print species amounts/concentrations for all orders
                if (params.allorders == 1):
                    conclist.append(linestart)
                else:
                    if (i == (len(data_block1)-1)):
                        conclist.append(linestart)
            
            # For interactions
            if (params.lbudget == 1):
                with open(os.path.join(tmpdir,"budget_" + str(icell)+"_"+str(timeperiod)+".pkl"),"rb") as fp:
                    for line in general_func.pickleLoader(fp):
                        data_block2.append(line)
                fp.close()
                for i in range(len(data_block2)-norders,len(data_block2)):
                    # Put riverid and cellnumber in front of the line and remove the time.
                    linestart = [rivernum[icell],rivmask[icell]]
                    linestart.extend(data_block2[i][1:])
                    budlist.append(linestart)    

            # For arguments
            if (params.largs == 1):
                with open(os.path.join(tmpdir,"arguments_" + str(icell)+"_"+str(timeperiod)+".pkl"),"rb") as fp:
                    for line in general_func.pickleLoader(fp):
                        data_block3.append(line)
                fp.close()
                for i in range(len(data_block3)-norders,len(data_block3)):
                    # Put riverid and cellnumber in front of the line and remove the time.
                    linestart = [rivernum[icell],rivmask[icell]]
                    linestart.extend(data_block3[i][1:])
                    argslist.append(linestart)         

        # Write all output of this riverbasin to the output file.
        stryearstart = general_func.stryear(yearstart)
        # Aquire lock
        file_locking.aquire_lock(params,lock,jobid,locktext="riverbasin"+stryearstart)

        fp = None
        if (params.allorders == 1):
            fp = open(os.path.join(params.outputdir,stryearstart+"_allorders.pkl"),"ab")
        else:
            fp = open(os.path.join(params.outputdir,stryearstart+".pkl"),"ab")
        for i in range(len(outputlist)):
            pickle.dump(outputlist[i],fp, -1)
        fp.close()

        fp = open(os.path.join(params.outputdir,"start"+stryearstart+".pkl"),"ab")
        for i in range(len(startoutputlist)):
            pickle.dump(startoutputlist[i],fp, -1)
        fp.close()

        fp = None
        if (params.allorders == 1):
            fp = open(os.path.join(params.outputdir,"conc_"+stryearstart+"_allorders.pkl"),"ab")
        else:
            fp = open(os.path.join(params.outputdir,"conc_"+stryearstart+".pkl"),"ab")

        for i in range(len(conclist)):
            pickle.dump(conclist[i],fp, -1)
        fp.close()

        if (params.lbudget == 1):
            fp = open(os.path.join(params.outputdir,"budget_"+stryearstart+".pkl"),"ab")
            for i in range(len(budlist)):
                pickle.dump(budlist[i],fp, -1)
            fp.close()

        if (params.largs == 1):
            fp = open(os.path.join(params.outputdir,"arguments_"+stryearstart+".pkl"),"ab")
            for i in range(len(argslist)):
                pickle.dump(argslist[i],fp, -1)
            fp.close()

        # Release lock
        file_locking.release_lock(params,lock,locktext="riverbasin"+stryearstart)


    # Now we can start with the dynamic calculations.
    starttime_calculate = time.time()
    
    print("Starting dynamic calculations...")
    timeperiod = -1
    lsteady = False
    yearstart = params.starttime
    yearend   = min(yearstart + params.outputtime,params.endtime)
    while not (yearstart > params.endtime):
        # Make a counter for each time period.
        timeperiod += 1
        lrestart = ((timeperiod+1)%iround(params.restartnumber) == 0) or\
                   (yearend == params.endtime)

        #lakes_dict = {}
        #outlakes_dict = {}
        # Loop over all cells:    
        for item in range(len(pointer1)):
            icell = pointer1[item].get_local_index()
            direction = lddmap.get_data(icell)
            next_cell = accuflux.goto_nextcell(icell,direction,lddmap.ncols,mask=rivmask)

            # Load is the time array of  all the sources            
            load = []
            for isrc in range(len(sources)):
                load.append(srcload[isrc][icell])       

            # Check whether cell is a lake/reservoir.
            lakeid = interpolate_list.stepwise(yearstart,lakeid_cells[icell],extrapol=1)[-1]
            if (lakeid > 0):
                outlakeid = interpolate_list.stepwise(yearstart,outlakeid_cells[icell],extrapol=1)[-1]
                endo_lake = interpolate_list.stepwise(yearstart,endo_lakes_cells[icell],extrapol=1)[-1]
                llake = True
                llakeout = (outlakeid  > 0.)
                lendolake = (endo_lake > 0.)
                #lakes_dict[icell] = iround(lakeid)
                #if (llakeout or lendolake):
                #    outlakes_dict[iround(lakeid)] = icell
            else:
                llake = False
                llakeout = False
                lendolake = False

            if (number_of_cpu > 1):
                #print "Submit nu ",icell
                try:
                    qq = lock_cell[next_cell]
                except KeyError:
                    lock_cell[next_cell] = mp.Lock()
            
                #if (prev_order != pointer1[item][0]):
                if (prev_order != pointer1[item].get_iorder()):
                    # We have to wait untill all processes are ready, before we can continue
                    for p in jobs:
                        p.join()
                    prev_order = pointer1[item].get_iorder()
                    jobs = []

                # Define the process with cellid "icell", return value p is process "handler"
                width = channel_width_grid.get_data(icell)
                slope = slopemap.get_data(icell)
                if (params.lstrahlergrids):
                    args_strahler_cell = args_strahler[icell]
                p=mp.Process(target= dgnm_cell.calculate_cell, args=(lsteady,lock_cell[next_cell],\
                                                     icell,params,species,proc,sources,\
                                                     tmpdir,next_cell,yearstart,yearend,timeperiod,\
                                                     runoff_cells[icell],load,temperature_cells[icell],\
                                                     globrad_cells[icell], discharge_cells[icell],\
                                                     volume_cells[icell],water_area_cells[icell],depth_cells[icell],\
                                                     width, slope, \
                                                     volume_fp_cells[icell],depth_fp_cells[icell],vel_cells[icell],\
                                                     low_veg_fr_cells[icell],high_veg_fr_cells[icell],llake,llakeout,lendolake,args_strahler_cell))
                
                # Add the process handler to the jobs list.
                jobs.append(p)
                # Start the process
                general_func.start_process(p,jobs,number_of_cpu)               
 
            else:
                # Do the job one at the time.
                width = channel_width_grid.get_data(icell)
                slope = slopemap.get_data(icell)
                if (params.lstrahlergrids):
                    args_strahler_cell = args_strahler[icell]
                dgnm_cell.calculate_cell(lsteady,None,icell,params,species,proc,sources,\
                                              tmpdir,\
                                              next_cell,yearstart,yearend,timeperiod,\
                                              runoff_cells[icell],load,temperature_cells[icell],\
                                              globrad_cells[icell], discharge_cells[icell],\
                                              volume_cells[icell],water_area_cells[icell],depth_cells[icell],\
                                              width, slope, \
                                              volume_fp_cells[icell],depth_fp_cells[icell],vel_cells[icell],\
                                              low_veg_fr_cells[icell],high_veg_fr_cells[icell],llake,llakeout,lendolake,args_strahler_cell=args_strahler_cell) 
            
        if (number_of_cpu > 1):            
            # Wait until all jobs are finished.
            for p in jobs:
                p.join()

        # Distribute the outflow of lake or reservoir to all the cells of the waterbody
        # Overwrite the information of the cells which are not the outflow of a waterbody.
        #for icell in lakes_dict:
        #    id = lakes_dict[icell]
        #    if (icell != outlakes_dict[id]):
        #        # This is a cell which must be overwritten with the outflow cell.
        #        outflowcell = outlakes_dict[id]
        #        my_sys.my_copyfile(os.path.join(tmpdir,str(outflowcell)+"_"+str(timeperiod)+".pkl"),\
        #                           os.path.join(tmpdir,str(icell)+"_"+str(timeperiod)+".pkl"))
        #        my_sys.my_copyfile(os.path.join(tmpdir,"conc_"+str(outflowcell)+"_"+str(timeperiod)+".pkl"),\
        #                           os.path.join(tmpdir,"conc_"+str(icell)+"_"+str(timeperiod)+".pkl"))
        #        if (params.lbudget == 1):
        #            my_sys.my_copyfile(os.path.join(tmpdir,"budget_"+str(outflowcell)+"_"+str(timeperiod)+".pkl"),\
        #                               os.path.join(tmpdir,"budget_"+str(icell)+"_"+str(timeperiod)+".pkl")) #LV 17-07-2017

        # Put cell information of the last time step into the world map.
        data_block3 = []
        data_block2 = []
        data_block1 = []
        data_block = []
        outputlist = []
        startoutputlist = []
        conclist = []
        budlist = []
        argslist = []
        for item in range(len(pointer1)):
            icell = pointer1[item].get_local_index()
            with open(os.path.join(tmpdir,str(icell)+"_"+str(timeperiod)+".pkl"),"rb") as fp:
                for line in general_func.pickleLoader(fp):
                    data_block.append(line)
            fp.close()
            # Store only the last output time of the cell.
            if (lrestart):
                for i in range(len(data_block)-norders,len(data_block)):
                    # Put riverid and cellnumber in front of the line and remove the time.
                    linestart = [rivernum[icell],rivmask[icell]]
                    linestart.extend(data_block[i][1:])
                    startoutputlist.append(linestart)

            for i in range(len(data_block)-norders,len(data_block)):
                # Put riverid and cellnumber in front of the line and remove the time.
                linestart = [rivernum[icell],rivmask[icell]]
                linestart.extend(data_block[i][1:])
                if (params.allorders == 1):
                    outputlist.append(linestart)
                else:
                    if (i == (len(data_block)-1)):
                        outputlist.append(linestart)

            with open(os.path.join(tmpdir,"conc_" + str(icell)+"_"+str(timeperiod)+".pkl"),"rb") as fp:
                for line in general_func.pickleLoader(fp):
                    data_block1.append(line)
            fp.close()

            for i in range(len(data_block1)-norders,len(data_block1)):
                # Put riverid and cellnumber in front of the line and remove the time.
                linestart = [rivernum[icell],rivmask[icell]]
                linestart.extend(data_block1[i][1:])
                if (params.allorders == 1):
                    conclist.append(linestart)
                else:
                    if (i == (len(data_block1)-1)):
                        conclist.append(linestart)

            if (params.lbudget == 1):
                with open(os.path.join(tmpdir,"budget_" + str(icell)+"_"+str(timeperiod)+".pkl"),"rb") as fp:
                    for line in general_func.pickleLoader(fp):
                        data_block2.append(line)
                fp.close()
                for i in range(len(data_block2)-norders,len(data_block2)):
                    # Put riverid and cellnumber in front of the line and remove the time.
                    linestart = [rivernum[icell],rivmask[icell]]
                    linestart.extend(data_block2[i][1:])
                    budlist.append(linestart)
                    #print(linestart)

            if (params.largs == 1):
                with open(os.path.join(tmpdir,"arguments_" + str(icell)+"_"+str(timeperiod)+".pkl"),"rb") as fp:
                    for line in general_func.pickleLoader(fp):
                        data_block3.append(line)
                fp.close()
                for i in range(len(data_block3)-norders,len(data_block3)):
                    # Put riverid and cellnumber in front of the line and remove the time.
                    linestart = [rivernum[icell],rivmask[icell]]
                    linestart.extend(data_block3[i][1:])
                    argslist.append(linestart)
					
        # Write all output of this riverbasin to the output file.
        # Aquire lock
        strtime = general_func.stryear(yearend)
        file_locking.aquire_lock(params,lock,jobid,locktext="riverbasin"+strtime)
        fp = None
        if (params.allorders == 1):
            fp = open(os.path.join(params.outputdir,strtime+"_allorders.pkl"),"ab")
        else:
            fp = open(os.path.join(params.outputdir,strtime+".pkl"),"ab")
        for i in range(len(outputlist)):
            pickle.dump(outputlist[i],fp, -1)
        fp.close()

        if (lrestart):
            fp = open(os.path.join(params.outputdir,"start"+strtime+".pkl"),"ab")
            for i in range(len(startoutputlist)):
                pickle.dump(startoutputlist[i],fp, -1)
            fp.close()

        fp = None
        if (params.allorders == 1):
            fp = open(os.path.join(params.outputdir,"conc_"+strtime+"_allorders.pkl"),"ab")
        else:
            fp = open(os.path.join(params.outputdir,"conc_"+strtime+".pkl"),"ab")
        for i in range(len(conclist)):
            pickle.dump(conclist[i],fp, -1)
        fp.close()

        if (params.lbudget == 1):
            fp = open(os.path.join(params.outputdir,"budget_"+strtime+".pkl"),"ab")
            for i in range(len(budlist)):
                pickle.dump(budlist[i],fp, -1)
            fp.close()

        if (params.largs == 1):
            fp = open(os.path.join(params.outputdir,"arguments_"+strtime+".pkl"),"ab")
            for i in range(len(argslist)):
                pickle.dump(argslist[i],fp, -1)
            fp.close()
		
        # Release lock
        file_locking.release_lock(params,lock,locktext="riverbasin"+strtime)

        print(" All processes are ready for river " + str(jobid) + " for time " + strtime+".")
        # Goto next periode.
        yearstart = yearend
        yearend   = min(yearstart + params.outputtime,params.endtime)
        if (abs(yearstart - yearend) < params.timestep):
            # This was the last timestep.
            break 

    # Remove the temporary directory with the output of this riverbasin
    my_sys.my_rmdir(tmpdir)

    endtime_calculate = time.time()
    print('Total time (in s) to calculate dynamic processes in riverbasin ' + str(iriv) + ':  ' + str(endtime_calculate-starttime_calculate))
                                                                                                     
if __name__ == '__main__':

    try:
        # Check command line arguments
        if (len(sys.argv) != 3):
            raise MyError("Not correct number of arguments on command line!",\
                          "Command line should be like: python dgnm_river.py general_objects.pkl rivers_1.pkl")

        # Get the command line arguments.
        parfile = sys.argv[1]
        numfile = sys.argv[2]

        # Check the input files
        if (not os.path.isfile(parfile)):
            raise MyError("Specified inputfile: " + parfile,\
                          "does not exist.")
        if (not os.path.isfile(numfile)):
            raise MyError("Specified inputfile: " + numfile,\
                          "does not exist.")

        # Get the general settings.   
        fp = open(parfile,"rb")
        params =  pickle.load(fp)
        species = pickle.load(fp)
        proc =    pickle.load(fp)
        sources = pickle.load(fp)
        fp.close()

        # Get river numbers and remove the input file.
        numlist = general_func.get_content_rm(numfile)

        # Run the model
        calculate(numlist,None,params,species,proc,sources)

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
        print("dgnm_river.py failed.")
        print(str(sys.exc_info()[0]))
        print(traceback.print_exc())
