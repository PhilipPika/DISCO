# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import general python modules
import os
import sys
import shutil
#import traceback
import pickle
import time as ts

import numpy as np
from netCDF4 import Dataset

# Import own general modules
import general_path
import ascraster
import directory
import manip
from error import MyError

# Import local modules
import general_func
#import analysis
import post_processing



def pkl_to_ascraster_allorders(filename,outputdir,unmask_raster,new_raster=None,label='',norders=6):
    # Convert output file in pkl format from rive_framework to ascraster format for each specie.
    # Open pkl file and read the header.
    #print('filename_outputconv=',filename)
    fp = open(filename,'rb')
    time,names = general_func.read_header(fp=fp)
    # Create for each specie an grid file.
    if (new_raster == None):
        duplicate = ascraster.duplicategrid(unmask_raster)
        lcoordinates_conversion = False
    else:
        duplicate = ascraster.duplicategrid(new_raster)
        lcoordinates_conversion = True

    grids=[]
    for iorder in range(norders):
        grids.append([ascraster.duplicategrid(duplicate)])
        for iname in range(1,len(names)):
            grids[-1].append(ascraster.duplicategrid(duplicate))

    # Now read the file line per line and put information into rasters
    try:
        order = 1
        while True:
            line = pickle.load(fp)
            # First element contains riverid, second the cell number in
            # unmasked environment and then for each specie a value.
            icell = line[1]

            # Get index in output raster.
            if (lcoordinates_conversion):
                #x,y = dummy.get_coordin_from_index(icell)
                x,y = unmask_raster.get_coordin_from_index(icell)
                icell = duplicate.get_index_from_coordin(x,y)

            # Set the values for each species
            for item in range(len(names)):
                grids[order-1][item].set_data(icell,line[item+2])

            order += 1
            if (order > norders):
                order = 1


    except EOFError:
        # End of file encountered.
        pass

    # Write all grids to disk
    for iorder in range(norders):
        for item in range(len(names)):
            strtime = str(round(time, 3))
            orderlabel = ''
            if (norders > 1):
                orderlabel = '_order'+str(iorder+1)
            grids[iorder][item].write_ascii_file(os.path.join(outputdir,os.path.join(strtime,label+names[item]+orderlabel+'.asc')))


def create_mask_from_pkl_file(filename, basinid):
    '''
    Returns a mask (numpy array of booleans)
    '''
    # Open pickle file
    fp = open(filename,'rb')
    time,names = general_func.read_header(fp=fp)

    # Create mask
    netcdf_mask = np.zeros((basinid.nrows, basinid.ncols), dtype=bool)
    netcdf_mask[:,:] = True

    # Now read the file line per line and get the cell numbers
    try:
        while True:
            line = pickle.load(fp)
            # First element contains riverid, second the cell number in unmasked environment
            icell =line[1]
            # Convert to ilat, ilon
            ilat, ilon = manip.calc_row_col_from_index(icell, basinid.ncols)

            netcdf_mask[ilat, ilon] = False

    except EOFError:
        # End of file encountered.
        pass

    return netcdf_mask


def init_ncdata(inputdir, ncdata, varname, basinid, unit='Mmol', long_name='None', fill_value=-999999999., zlib=True,\
                institution='Dept. of Earth Sciences - Geochemistry, Utrecht University',\
                description='Conversion of DGNM output with output_conversion.py'):
    '''
    Initializes NetCDF file and fills metadata and creates dimensions of NETCDF data
    '''

    ncdata.institution = institution
    ncdata.description = description

    #### Create dimensions: time is unlimited, others are fixed
    ncdata.createDimension('time', None)
    ncdata.createDimension('lat', basinid.nrows)
    ncdata.createDimension('lon', basinid.ncols)

    #### Create variables
    date_time = ncdata.createVariable('time','f4',('time',))
    date_time.standard_name = 'time'
    date_time.units = 'hours since 1900-01-01 00:00:00'
    date_time.calendar = 'standard'

    lat = ncdata.createVariable('lat','f4',('lat',))
    lat.standard_name = 'latitude'
    lat.long_name = 'Latitude of cell centers'
    lat.units = 'degrees north'
    lat[:] = manip.initialize_lat(basinid.nrows, basinid.yllcorner, basinid.cellsize)

    lon = ncdata.createVariable('lon','f4',('lon',))
    lon.standard_name = 'longitude'
    lon.long_name = 'Longitude of cell centers'
    lon.units = 'degrees east'
    lon[:] = manip.initialize_lon(basinid.ncols, basinid.xllcorner, basinid.cellsize)

    var = ncdata.createVariable(varname,'f4',('time','lat','lon',) ,fill_value=fill_value,zlib=zlib, complevel=4)
    var.standard_name = varname
    if (varname.startswith('conc_')):
        var.long_name = 'Concentration of ' + varname[5:]
        var.units = 'mg per liter'
    elif varname=='pH':
        var.long_name = 'Acidity'
        var.units = '-log of [H+]'
    elif varname=='pCO2':
        var.long_name = 'CO2 pressure'
        var.units = 'parts per million'
    elif unit=='Mmol':
        var.long_name = 'Amount of ' + varname[5:]
        var.unit = unit
    elif 'loadIN' in varname:
        v_name = varname[:varname.find('loadIN')]
        vdum = str()
        for j in [i for i, c in enumerate(v_name) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = "imported "+vdum+" from upstream"
    elif 'loadOUT' in varname:
        v_name = varname[:varname.find('loadOUT')]
        vdum = str()
        for j in [i for i, c in enumerate(v_name) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = "exported "+vdum+" to downstream"
    elif '_dvoldt' in varname:
        vdum = str()
        for j in [i for i, c in enumerate(varname) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = 'change of '+vdum+' due to volume change'
    elif 'atm_exch' in varname:
        v_name = varname[:varname.find('atm_exch')]
        vdum = str()
        for j in [i for i, c in enumerate(v_name) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = "exchange of "+vdum+" with atmosphere"
    elif 'doc_excretion' in varname:
        vdum = str()
        for j in [i for i, c in enumerate(varname) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = "excretion of DOC by "+vdum
    elif 'oxidation' in varname:
        spec = varname[:-3].replace("_oxidation_order", "")
        var.long_name = "oxidation of "+spec.upper()
    elif 'grazing' in varname:
        vdum = str()
        for j in [i for i, c in enumerate(varname) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = "grazing by "+vdum
    elif 'mortality' in varname:
        vdum = str()
        for j in [i for i, c in enumerate(varname) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = "mortality of "+vdum
    elif 'prim_prod' in varname:
        vdum = str()
        for j in [i for i, c in enumerate(varname) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = "primary production by "+vdum
    elif 'respiration' in varname:
        vdum = str()
        for j in [i for i, c in enumerate(varname) if c.isupper()]:
            vdum+=str(varname[j])
        var.long_name = "respiration by "+vdum
    else:
        var.long_name = long_name
        var.unit = unit

def pkl2netcdf_allorders(filename,spec_data,netcdf_mask,basinid,fill_value=-999999999.,label="",norders=6):
    '''
    Convert output file in pkl format from rive_framework to NETCDF format for each specie.
    Open pkl file and read the header.
    '''

    fp = open(filename,'rb')
    time,names = general_func.read_header(fp=fp)

    # Create empty masked grid for each species
    netcdf_grids=[]
    for iorder in range(norders):
        netcdf_grids.append([])
        #for iname in range(1,len(names)):
        for iname in range(len(names)):
            netcdf_grids[-1].append(np.ma.array(np.zeros((basinid.nrows,basinid.ncols)), mask=netcdf_mask, fill_value=fill_value).unshare_mask())

    # Now read the file line per line and put information into rasters
    try:
        order = 1
        while True:
            line = pickle.load(fp)
            # First element contains riverid, second the cell number in unmasked environment and then for each specie a value.
            icell =line[1]
            ilat, ilon = manip.calc_row_col_from_index(icell, basinid.ncols)
            # Set the values for each species
            for item in range(len(names)):
                netcdf_grids[order-1][item][ilat,ilon] = line[item+2]

            order += 1
            if (order > norders):
                order = 1
        fp.close()
    except(EOFError):
        # End of file encountered.
        pass

    time_unit = 'hours since 1900-01-01 00:00:00'
    t0 = manip.convert_year2numdate(float(time), time_unit)
    for iorder in range(norders):
        for item in range(len(names)):
            orderlabel = ''
            if (norders > 1):
                orderlabel = '_order'+str(iorder+1)
            varname = names[item]

            if (label == "conc_"):
                varname = label+varname
            manip.add_grid_time(spec_data[iorder][item],varname,netcdf_grids[iorder][item],t0)





def Convert_concentration(outputfiles, filename, basinid, params, conc_outputfiles):
    #### Convert species concentration
    print("### Conversion of the species concentration... ###")
    starttime_Conc = ts.time()
    if len(outputfiles) > 0:
        filename = outputfiles[0][0]
        # Read header of one output file to get species names
        time,names = general_func.read_header(filename=filename)

        # Create mask for NETCDF data
        netcdf_mask = create_mask_from_pkl_file(filename, basinid)

        # Use speciesnames to create NETCDF datasets in output directory
        conc_spec_data = []
        for name in names: # Create NETCDF output files
            states_dir = directory.ensure(os.path.join(params.outputdir, '..', 'STATES'))
            ncfile = os.path.join(states_dir,'conc_'+name+'.nc')
            conc_spec_data.append(Dataset(ncfile,'w'))
            # Initialize each of the NETCDF files (metadata, dimensions...)
            init_ncdata(inputdir,conc_spec_data[-1], 'conc_'+name, basinid)

        # Fill in the NETCDF concentration output files
        for filename,time in conc_outputfiles:
            # Conversion of pkl file to ascraster files for each speciename
            pkl2netcdf_allorders(filename,[conc_spec_data],netcdf_mask,basinid,label='conc_',norders=1)

        for item in range(len(names)):
            conc_spec_data[item].sync()
            conc_spec_data[item].close()
    endtime_Conc = ts.time()
    print('...finished in (in s):  ' + str(endtime_Conc-starttime_Conc))



def Convert_concentration_sub(conc_outputfiles_allorders, filename, basinid, params, outputfiles_allorders, norders):
    #### Convert species concentration in subgrid orders
    print("### Conversion of the species concentration in subgrid orders ... ###")
    starttime_ConcSub = ts.time()
    if len(outputfiles_allorders) > 0:
        filename = outputfiles_allorders[0][0]
        # Read header of one output file to get species names
        time,names = general_func.read_header(filename=filename)

        # Create mask for NETCDF data
        netcdf_mask = create_mask_from_pkl_file(filename, basinid)


        # Calculates the concentration of the subgrid orders
        conc_spec_data_allorders = []
        for iorder in range(norders):
            conc_spec_data_allorders.append([])
            orderlabel = '_order'+str(iorder+1)

            for name in names: # Create concentration files for subgrid orders
                subgrid_states_dir = directory.ensure(os.path.join(params.outputdir, '..', 'STATES', 'subgrid'))
                ncfile = os.path.join(subgrid_states_dir,'conc_'+name+orderlabel+'.nc')
                conc_spec_data_allorders[-1].append(Dataset(ncfile,'w'))
                # Initialize each of the NETCDF files (metadata, dimensions...)
                init_ncdata(inputdir, conc_spec_data_allorders[-1][-1], 'conc_'+name, basinid)

        for filename,time in conc_outputfiles_allorders:
            pkl2netcdf_allorders(filename,conc_spec_data_allorders,netcdf_mask,basinid,label='conc_',norders=norders)

        for iorder in range(norders):
            for item in range(len(names)):
                conc_spec_data_allorders[iorder][item].sync()
                conc_spec_data_allorders[iorder][item].close()

        for fn in directory.get_files_with_str(subgrid_states_dir, "*order6*"):
            shutil.copyfile(fn, os.path.join(subgrid_states_dir, '..', os.path.basename(fn).replace("_order6", "")))
    endtime_ConcSub = ts.time()
    print('...finished in (in s):  ' + str(endtime_ConcSub-starttime_ConcSub))



def Convert_env_hydro(argumentfiles, filename, basinid, params, norders):
    #### Convert environmental parameters and hydrology
    print("### Conversion of the environmental parameters and hydrology ... ###")
    starttime_Env = ts.time()
    if (len(argumentfiles)>0):
        filename = argumentfiles[0][0]
        # Read header of one output file to get argument names
        time,names = general_func.read_header(filename=filename)

        # Create mask for NETCDF data
        netcdf_mask = create_mask_from_pkl_file(filename, basinid)

        argument_data = []

        for iorder in range(norders):
            argument_data.append([])

            orderlabel = '_order'+str(iorder+1)
            for name in names: # Create NETCDF output files
                subgrid_arg_dir = directory.ensure(\
                  os.path.join(params.outputdir, '..', 'STREAM_ENV_CONDITIONS', 'subgrid'))
                ncfile = os.path.join(subgrid_arg_dir,name+orderlabel+'.nc')
                argument_data[-1].append(Dataset(ncfile,'w'))
                # Initialize each of the NETCDF files (metadata, dimensions...)
                u = 'None'
                long_name = 'None'
                ln = long_name
                if name=='vol' or 'vol_' in name:
                    u='[km3]'
                    ln='volume of water'
                elif name=='depth' or 'depth_' in name:
                    u='[m]'
                    ln='depth of waterbody'
                elif name=='discharge' or 'discharge_' in name:
                    u='[km3/yr]'
                    ln='discharge of waterbody'
                elif name=='dvoldt' or 'dvoldt_' in name:
                    u='[km3/yr]'
                    ln='volume change of waterbody'
                elif name=='flow_velocity' or 'flow_velocity_' in name:
                    u='[km/yr]'
                    ln='flow velocity of waterbody'
                elif name=='globrad' or 'globrad_' in name:
                    u='[W/m2]'
                    ln='global radiation upon waterbody'
                elif name=='length' or 'length_' in name:
                    u='[m]'
                    ln='length of waterbody'
                elif name=='residence_time' or 'residence_time_' in name:
                    u='[yr]'
                    ln='residence time of the waterbody'
                elif name=='slope' or 'slope_' in name:
                    u='[unknown unit]'
                    ln='slope of the waterbody'
                elif name=='temperature' or 'temperature_' in name:
                    u='[K]'
                    ln='temperature of the waterbody'
                elif name=='width' or 'width_' in name:
                    u='[m]'
                    ln='width of the waterbody'
                elif name=='windspeed' or 'windspeed_' in name:
                    u='[m/s]'
                    ln='windspeed 10 meter above waterbody'
                elif name=='area' or 'area_' in name:
                    u='[km2]'
                    ln='area of the waterbody'
                elif name=='icell' or 'icell_' in name:
                    u='[#]'
                    ln='cell id within mask'
                init_ncdata(inputdir, argument_data[-1][-1], name, basinid, unit=u, long_name=ln)

        for filename,time in argumentfiles:
            pkl2netcdf_allorders(filename,argument_data,netcdf_mask,basinid, norders=norders)

        for iorder in range(norders):
            for item in range(len(names)):
                argument_data[iorder][item].sync()
                argument_data[iorder][item].close()

        for fn in directory.get_files_with_str(subgrid_arg_dir, "*order6*"):
            shutil.copyfile(fn, os.path.join(subgrid_arg_dir, '..', os.path.basename(fn).replace("_order6", "")))
    endtime_Env = ts.time()
    print('...finished in (in s):  ' + str(endtime_Env-starttime_Env))



def Convert_budget(budgetfiles, filename, basinid, params, norders):
    #### Convert of budgets
    print("### Conversion of budgets ... ###")
    starttime_Spec = ts.time()


    if (len(budgetfiles)>0):
        filename = budgetfiles[0][0]
        # Read header of one output file to get argument names
        time,names = general_func.read_header(filename=filename)
        # Create mask for NETCDF data
        netcdf_mask = create_mask_from_pkl_file(filename, basinid)

        # Create converted NETCDF output files for each species
        # Use speciesnames to create NETCDF datasets in output directory
        budget_data = []

        # starttime_iorder = ts.time()
        for iorder in range(norders):

            # print('Budget Order: ', iorder)
            budget_data.append([])
            orderlabel = '_order'+str(iorder+1)

            for name in names: # Create NETCDF output files
                # starttime_name = ts.time()
                budget_dir = directory.ensure(os.path.join(params.outputdir, '..', 'BUDGET', 'subgrid'))
                ncfile = os.path.join(budget_dir,name+orderlabel+'.nc')

                budget_data[-1].append(Dataset(ncfile,'w'))
                # Initialize each of the NETCDF files (metadata, dimensions...)
                init_ncdata(inputdir, budget_data[-1][-1], name, basinid)
                # endtime_name = ts.time()
                # print('name lasted (in s):  ' + str(endtime_name-starttime_name))

        starttime_iorder = ts.time()
        # LV 12-07-2017
        for filename,time in budgetfiles:
            # print('Budgetfile Name: ', filename)
            pkl2netcdf_allorders(filename,budget_data,netcdf_mask,basinid,label="budget_", norders=norders)
        endtime_order = ts.time()
        print('iorder lasted (in s):  ' + str(endtime_order-starttime_iorder))


        # Parallelizing with Pool.starmap()
        # import multiprocessing as mp
        # pool = mp.Pool(mp.cpu_count())
        # [pool.apply(pkl2netcdf_allorders, \
        #     args = (budgetfiles,budget_data,netcdf_mask,basinid,"budget_", norders))]
        # pool.close()


        for iorder in range(norders):
            for item in range(len(names)):
                budget_data[iorder][item].sync()
                budget_data[iorder][item].close()

        for fn in directory.get_files_with_str(budget_dir, "*order6*"):
            # print('copy file: ', fn)
            shutil.copyfile(fn, os.path.join(budget_dir, '..', os.path.basename(fn).replace("_order6", "")))



        endtime_Spec = ts.time()
        print('...finished in (in s):  ' + str(endtime_Spec-starttime_Spec))





def convert_output(inputdir,outformat,norders=6):
    #### Initial setup - file containers
    print("Entering output_conversion >> convert_output()")
    '''
    Read output of rive model and create output grids in wanted format (ASCII, NETCDF or both).
    By default norders is set to Strahler order = 6
    '''

    # Check if the format conversion asked by the user is possible
    if not((outformat == 'ASCII') or (outformat == 'NETCDF') or (outformat == 'ALL')):
        print('Unknown conversion format provided to convert_output function: ' + outformat)
        print('Only ASCII, NETCDF or ALL are accessible')
        raise MyError('Run again with one of these 3 formats')

    if inputdir.startswith('..'):
        inputdir =  os.path.join(os.getcwd(), inputdir)

    # Get the information of the mask to an output file to be able to read the output files.
    pklFile = open(os.path.join(inputdir,'mask.pkl'),'rb')
    try:
        basinid = pickle.load(pklFile)
        params =  pickle.load(pklFile)
    except KeyError:
        raise MyError("Something wrong with reading the meta information of the output directory.")
    pklFile.close()

    if params.lfloodplains==1:
        norders+=1

    # Make the output directory
    #if (not os.path.isdir(outputdir)):
    #    os.makedirs(outputdir)

    # Get all files from the input directory
    listdir = os.listdir(inputdir)

    outputfiles = []
    conc_outputfiles = []
    outputfiles_allorders = []
    conc_outputfiles_allorders = []
    argumentfiles = []
    budgetfiles = []

    for filename in listdir:
        fullname = os.path.join(inputdir,filename)
        if os.path.isfile(fullname):
            # Split filename into extension and basename
            basename,extension = os.path.splitext(filename)

            # Copy the log file of the simulation run to the new output directory
            #if (extension.upper() == '.LOG'):
            #    shutil.copy(fullname,os.path.join()

            # Select all the files which have the extension ".pkl"
            if extension.upper() == '.PKL':
                # Try whether the basename is a float (time) or starts by "conc_"
                try:
                    # addition LV 12-07-2017: option to print outputs for all orders
                    if basename.endswith("_allorders"):
                        if basename.startswith("conc_"):
                            time = float(basename[5:-10])
                            conc_outputfiles_allorders.append([fullname,time])
                        else:
                            time = float(basename[:-10])
                            outputfiles_allorders.append([fullname,time])
                    else:
                        if (not basename.startswith("conc_") and not \
                            basename.startswith("arguments_") and not \
                              "budget_" in basename):
                            time = float(basename)
                            # This is a file which we want to convert.
                            outputfiles.append([fullname,time])
                        elif basename.startswith("conc_"):
                            time = float(basename[5:])
                            # This is a file which we want to convert.
                            conc_outputfiles.append([fullname,time])
                        elif basename.startswith("arguments_"):
                            time = float(basename[10:])
                            # This is a file which we want to convert.
                            argumentfiles.append([fullname,time])
                        elif basename.startswith("budget_"):
                            time = float(basename[7:])
                            # This is a file which we want to convert.
                            budgetfiles.append([fullname,time])
                except ValueError:
                    pass
    # Sorting the list of input data
    outputfiles.sort()
    conc_outputfiles.sort()
    outputfiles_allorders.sort()
    conc_outputfiles_allorders.sort()
    argumentfiles.sort()
    budgetfiles.sort()

        # Conversion to ASCII - var 'outputdir' not undefined - potential error?
    if ((outformat == 'ASCII') or (outformat == 'ALL')):

        # Create a dummy grid which is not masked. This is used to convert coordinates.
        dummy = ascraster.Asciigrid(ncols=basinid.ncols,nrows=basinid.nrows,\
                                    xllcorner=basinid.xllcorner,yllcorner=basinid.yllcorner,\
                                    cellsize=basinid.cellsize,nodata_value=-9999.)

        # Use output dates to create subdirectories in output directory
        for item in range(len(outputfiles)):
            strtime = str(round(outputfiles[item][-1], 3))
            if (not os.path.isdir(os.path.join(params.outputdir,strtime))):
                os.makedirs(os.path.join(params.outputdir,strtime))

        # Conversion of pkl file to ascraster file for each output date and each species
        for filename,time in outputfiles:
            pkl_to_ascraster_allorders(filename,params.outputdir,dummy,norders=1)

        # Conversion of pkl file to ascreaster of concentration files
        for filename,time in conc_outputfiles:
            pkl_to_ascraster_allorders(filename,params.outputdir,dummy,\
                                       label='conc_',norders=1)

        # LV 12-07-2017
        for item in range(len(outputfiles_allorders)):
            strtime = str(round(outputfiles_allorders[item][-1], 3))
            if (not os.path.isdir(os.path.join(params.outputdir,strtime))):
                os.makedirs(os.path.join(params.outputdir,strtime))

        for filename,time in outputfiles_allorders:
            pkl_to_ascraster_allorders(filename,params.outputdir,dummy)

        for filename,time in conc_outputfiles_allorders:
            pkl_to_ascraster_allorders(filename,params.outputdir,dummy,label='conc_')

        for filename,time in argumentfiles:
            pkl_to_ascraster_allorders(filename,params.outputdir,dummy,label='arguments_')





    if ((outformat == 'NETCDF') or (outformat == 'ALL')):

        if (len(outputfiles) == 0) and (len(outputfiles_allorders) == 0)\
            and (len(argumentfiles) ==0) and (len(budgetfiles) == 0):
            raise MyError("There is nothing to do!")




    #### Convert output - calling functions
    Convert_concentration(outputfiles, filename, basinid, params, conc_outputfiles)

    Convert_concentration_sub(conc_outputfiles_allorders, filename, basinid, \
                              params, outputfiles_allorders, norders)

    Convert_env_hydro(argumentfiles, filename, basinid, params, norders)

    Convert_budget(budgetfiles, filename, basinid, params, norders)

    #### Calcaulte diagnostic variables pH and pCO2
    post_processing.calc_ph_pco2(params, outformat) # PP This script converts pH and pCO2 for the given data






if __name__ == "__main__":
    # Set the general path for the own python modules

    if (len(sys.argv) < 2):
        print("Not enough arguments on command line!")
        raise MyError("Command line should be: python output_conversion.py inputdir outformat")
    inputdir = sys.argv[1]
    outformat = sys.argv[2]

    print("OUTPUT CONVERSION OF " , inputdir, " TO ", outformat)

    starttime_main = ts.time()

    convert_output(inputdir,outformat)

    endtime_main = ts.time()
    print('Total OUTPUT CONVERSION lasted (in s):  ' + str(endtime_main-starttime_main))