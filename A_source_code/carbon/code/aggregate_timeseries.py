# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
import copy
import os
import sys


from netCDF4 import Dataset
#import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
#import imageio
#import csv
import pandas as pd

import general_path

import read_parameter

import ascraster
#import cmd_options_dgnm
import directory
#import get_dynamic_gridinfo
import general_startup
#import interpolate_list
import make_index_species
#import make_outlakes_dict
import make_mask
# import output_conversion
# import pickle
#import pointer_class
import reactions
#import specie
import define_subgrid_streamorder

import manip
import time


global_colors = cm.Set1(np.linspace(0,1,8))

def do(params):
    starttime_main_all = time.time()
    ## make all tables
    starttime_main = time.time()
    all_inputs_to_table(params)
    endtime_main = time.time()
    print('all_inputs_to_table finished (in s):  ' + str(endtime_main-starttime_main))

    starttime_main = time.time()
    all_atm_exch_to_table(params)
    endtime_main = time.time()
    print('all_atm_exch_to_table finished (in s):  ' + str(endtime_main-starttime_main))

    starttime_main = time.time()
    all_exports_to_table(params)
    endtime_main = time.time()
    print('all_exports_to_table finished (in s):  ' + str(endtime_main-starttime_main))

    starttime_main = time.time()
    all_budget_to_table(params)
    endtime_main = time.time()
    print('all_budget_to_table finished (in s):  ' + str(endtime_main-starttime_main))

    starttime_main = time.time()
    all_fluxes_to_table(params)
    endtime_main = time.time()
    print('all_fluxes_to_table finished (in s):  ' + str(endtime_main-starttime_main))

    starttime_main = time.time()
    conv_all_tables_to_Tg(params)
    endtime_main = time.time()
    print('conv_all_tables_to_Tg finished (in s):  ' + str(endtime_main-starttime_main))

    print('aggregate timeseries.py finished (in s):  ' + str(endtime_main-starttime_main_all))

    # all_stream_env_conditions_to_table(params) # function exists but unused

def get_river_name(params):
    if (not 'country' in params.file_mask) and (params.mask_bool_operator=='EQ'):
        if params.maskid==1:
            rivername = 'Amazon'
        elif params.maskid==2:
            rivername = 'Congo'
        elif params.maskid==3:
            rivername = 'Caspian Sea'
        elif params.maskid==4:
            rivername = 'Mississippi'
        elif params.maskid==5:
            rivername = 'Nile'
        elif params.maskid==7:
            rivername = 'Yenisei'
        elif params.maskid==8:
            rivername = 'Ob'
        elif params.maskid==9:
            rivername = 'Lena'
        elif params.maskid==10:
            rivername = 'Yangtze'
        elif params.maskid==11:
            rivername = 'Amur'
        elif params.maskid==13:
            rivername = 'Mackenzie'
        elif params.maskid==28:
            rivername = 'Yukon'
        elif params.maskid==38:
            rivername = 'Kolyma'
        elif params.maskid==55:
            rivername = 'Indigirka'
        elif params.maskid==99:
            rivername = "Rhine"
        elif params.maskid==220:
            rivername = 'Peel'
        elif params.maskid==231:
            rivername = 'Seine'
        elif params.maskid==433:
            rivername = 'Meuse'
        elif params.maskid==4634:
            rivername = 'Shanghai'
        else:
            rivername = str(params.maskid)
    elif params.maskid==0 and params.mask_bool_operator=='GT':
        rivername = 'Global'
    elif 'country' in params.file_mask:
        if params.maskid==752:
            rivername = "Sweden"
    else:
        rivername = 'other'
    return rivername

def geographical_reach_npmask(params):
    mincol = 1e6
    maxcol = 0
    minrow = 1e6
    maxrow = 0
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')

    for row in range(dum_asc.nrows):
        for col in range(dum_asc.ncols):
            if mask[row,col]==False:
                if col > maxcol:
                    maxcol = col
                if col < mincol:
                    mincol = col
                if row > maxrow:
                    maxrow = row
                if row < minrow:
                    minrow = row
    minlon, maxlat = dum_asc.get_coordin_from_row_col(minrow,mincol)
    maxlon, minlat = dum_asc.get_coordin_from_row_col(maxrow,maxcol)

    d = dict()
    d['minlon']=minlon
    d['maxlon']=maxlon
    d['minlat']=minlat
    d['maxlat']=maxlat
    return d

def debugprint(params, print_statement):
    if params.ldebug: print(print_statement)


def make_3d_mask(params, species, folder, mask_kind):
    # changes here: folder
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")
    #invariant
    dummy_nc = Dataset(proclist[-1], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[-1]))[0][:-7]

    modeldat_startindex, modeldat_endindex, wbody_startindex, wbody_endindex, waterbodyid = \
    make_time_indices(params, dummy_nc)


    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    # chnages here
    if mask_kind == 'climate.asc':
        mask_fn = os.path.join(params.water_inputdir, 'climate.asc')
        mask_2d_dum = make_mask.do(mask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')

        # invariant: mask
        mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc,logical = params.mask_bool_operator, mask_type='np_grid')
        mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
        mask_2d[:, :] = True
        mask_2d[np.where(np.logical_and(mask_2d_dum[:, :] == False,mask_2d_dum[:, :] == False))] = False

    elif mask_kind=='rivermouth.asc': # rivermouth
        mask_fn = os.path.join(params.water_inputdir, "rivermouth.asc")
        dum_mask = ascraster.create_mask(
        mask_fn, params.maskid, logical=params.mask_bool_operator, numtype=int)
        if len(dum_mask) > 0:
            mask_fn = os.path.join(params.water_inputdir, "rivermouth.asc")
        else:  # endoreic basin
            mask_fn = os.path.join(params.water_inputdir,"rivermouths_exporting_to_endoreic_lakes.asc")
        mask_2d = make_mask.do(mask_fn, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    else:
        print('Mask should be either rivermouth.asc OR climate.asc')

    #invariant
    mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)
    return mask_3d, mask_2d, dummy_name, dummy_nc, dum_asc,  modeldat_startindex, modeldat_endindex, wbody_startindex, wbody_endindex, waterbodyid


def make_time_indices(params, dummy_nc):
    if params.outputtime < 1:
        print("MONTHLY DATA - from waterbodyid_2399months")
        debugprint(params,'Time step')
        debugprint(params,params.outputtime)
        debugprint(params,"find index for outputtime ")
        debugprint(params,'size dummy nc')
        debugprint(params,'print dummy nc')
        debugprint(params,dummy_nc)
        debugprint(params,'print dummy nc TIME')
        debugprint(params,dummy_nc['time'][:])

        debugprint(params,os.path.join(params.water_inputdir))
       #This parameter 'waterbodyid' needs to have the same hourly dt (time step) as the DISC model output.
       # In my case hours/month, while chosing days/month, as with the new hydrology input, won't work.

        # waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_1970_2010.nc"), 'r')
        waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_2399months.nc"), 'r')

        debugprint(params,'Print waterbody id time vector')
        debugprint(params,waterbodyid['time'][:])


        modelrun_dummy_start = max(dummy_nc['time'][0],0)
        debugprint(params,'modelrun_dummy_start')
        debugprint(params,modelrun_dummy_start)


        modelrun_dummy_end = dummy_nc['time'][-1]
        debugprint(params,'modelrun_dummy_end')
        debugprint(params,modelrun_dummy_end)

        debugprint(params,'TESTing the time step overlap')
        debugprint(params,np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0][0])
        wbody_startindex = np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0][0]
        wbody_startindex=max(wbody_startindex, 0)
        debugprint(params,'All data start index')
        debugprint(params,wbody_startindex)
        debugprint(params,"wbody_startindex from params,np.where(waterbodyid['time'][:] >= modelrun_dummy_start)")
        debugprint(params, np.where(waterbodyid['time'][:] >= modelrun_dummy_start))
         # One [0] shows the array in position 1 aka (:,:)
        debugprint(params,"wbody_startindex from np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0]")
        debugprint(params, np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0])
         # two [0] show the first nr of first dim = [1,1]
        debugprint(params,"wbody_startindex from params,np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0][0]")
        debugprint(params, np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0][0])
        debugprint(params,"wbody_startindex from params,np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0][-1]+1")
        debugprint(params, np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0][-1]+1)

        wbody_endindex = np.where(waterbodyid['time'][:] <= modelrun_dummy_end)[0][-1]+1

        debugprint(params,'All data end index')
        debugprint(params,wbody_endindex)

        debugprint(params,'waterbodyid time steps: -1, 0, 1')
        debugprint(params,waterbodyid['time'][-1])
        debugprint(params,waterbodyid['time'][0])
        debugprint(params,waterbodyid['time'][1])

        modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyid['time'][wbody_startindex])[0][0]
        debugprint(params,'Overlap In-/output start index')
        debugprint(params,modeldat_startindex)



        modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyid['time'][wbody_endindex])[0][-1]+2
        debugprint(params,'Overlap In-/output end index')
        debugprint(params,modeldat_endindex)

        if modeldat_endindex==len(dummy_nc['time'][:]):
            wbody_endindex +=1
    else:
        print("MONTHLY DATA - from waterbodyid_201yrs")
        waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_201yrs.nc"), 'r')
        modelrun_dummy_start = max(dummy_nc['time'][0],0)
        modelrun_dummy_end = dummy_nc['time'][-1]
        wbody_startindex = np.where(waterbodyid['time'][:] >= modelrun_dummy_start)[0][0]
        wbody_startindex=max(wbody_startindex, 0)
        wbody_endindex = np.where(waterbodyid['time'][:] <= modelrun_dummy_end)[0][-1]
        modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyid['time'][wbody_startindex])[0][0]
        modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyid['time'][wbody_endindex])[0][-1]+1
        if modeldat_endindex==len(dummy_nc['time'][:]):
            wbody_endindex +=1
    return modeldat_startindex, modeldat_endindex, wbody_startindex, wbody_endindex, waterbodyid


def dict_to_csv(filename, mydict):
    pd.DataFrame(mydict).T.reset_index().to_csv(filename, header=False, index=False)

def csv_to_dict(filename):
    dummy_dict = pd.read_csv(filename, index_col=0, header=None).to_dict("split")
    mydict = dict(zip(dummy_dict['index'], dummy_dict['data']))
    return mydict

def conv_all_tables_to_Tg(params):
    table_dir = os.path.join(params.outputdir, "..", "ANALYSIS", "tables")
    table_list = directory.get_files_with_str(table_dir, "*.csv", exclude=['validation', 'Tg'])
    conv = 12*(1/params.outputtime)*1e-6
    for table_filename in table_list:
      flux_dict = csv_to_dict(table_filename)
      new_flux_dict = {}
      for key, array in flux_dict.items():
        if not key=='time':
          new_flux_dict[key] = (np.array(array)*conv).tolist()
        else:
          new_flux_dict[key] = array
      new_table_filename = table_filename[:-4]+'_Tg'+table_filename[-4:]
      dict_to_csv(new_table_filename, new_flux_dict)

def all_inputs_to_dict(params,add_color=False):

    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")

    mask_3d, mask_2d, dummy_name, dummy_nc, dum_asc,  modeldat_startindex, \
        modeldat_endindex, wbody_startindex, wbody_endindex, waterbodyid= \
        make_3d_mask(params, species, folder, 'climate.asc')

    src_series = dict()
    src_series["time"] = manip.convert_numdate2year(dummy_nc['time'][modeldat_startindex:modeldat_endindex], dummy_nc['time'].units)

    tot_in = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))

#    for specie in species:
#      src_files = sorted(directory.get_files_with_str(folder, "atmospheric_exchange_"+specie.get_name().upper()+'*.nc'))

    ## this is still a bit too specific for carbon; to be formulated more generic
    for source in sources:
      src_nc = Dataset(os.path.join(params.load_inputdir, source.get_val('name')+'.nc'), 'r')
      debugprint(params,'Current src_nc read:')
      debugprint(params,src_nc)
      for attrib in source.get_attrib():
          if 'fr_' in attrib:
            if (not 'TSS' in attrib) and (not 'PIM' in attrib):
              fraction = getattr(source, attrib)
            elif ('tss' in source.get_val('name').lower() and ('TSS' in attrib)) or ('pim' in source.get_val('name').lower() and ('PIM' in attrib)):
              fraction = 1

      debugprint(params,source.get_val('name'))
      A=src_nc[source.get_val('name')]
      debugprint(params,'src_nc shape pre modification')
      debugprint(params,A.shape)
      debugprint(params, wbody_startindex)
      debugprint(params, wbody_endindex)
      src_grid = src_nc[source.get_val('name')][wbody_startindex:wbody_endindex,:,:]*params.outputtime*fraction

      debugprint(params,'src_nc shape post modification')
      debugprint(params,src_grid.shape)
      debugprint(params,'Mask size')
      debugprint(params,mask_3d.shape)

      # src_series[source.get_val('name')] = np.nansum(np.nansum(np.ma.array(src_grid, mask=mask_3d),axis=2),axis=1).tolist()
      debugprint(params,'FINISHED SOURCE MULTIPLICATION')


      for attrib in source.get_attrib():
        if 'fr_' in attrib:
          if (not 'TSS' in attrib) and (not 'ALK' in attrib) and (not 'PIM' in attrib):
            tot_in = np.add(src_grid, tot_in)
      for attrib in source.get_attrib():
        if 'fr_' in attrib:
          src_grid = src_nc[source.get_val('name')][wbody_startindex:wbody_endindex,:,:]*params.outputtime*fraction

          if (not 'TSS' in attrib) and (not 'PIM' in attrib) :
            source_string = attrib.replace('fr_', '')+'_srcloadIN'
          elif 'tss' in source.get_val('name').lower():
            source_string = attrib.replace('fr_', '')+'_srcloadIN'
          elif ('tss' not in source.get_val('name').lower()) and ('TSS' in attrib):
            source_string = 'TSS_srcloadIN'
            src_grid*=source.get_val(attrib)
          elif ('pim' not in source.get_val('name').lower()) and ('PIM' in attrib):
            source_string = 'PIM_srcloadIN'
            src_grid*=source.get_val(attrib)
          if source_string in list(src_series.keys()):
            src_series[source_string] = np.add(np.array(src_series[source_string]), np.nansum(np.nansum(np.ma.array(src_grid, mask=mask_3d),axis=2),axis=1)).tolist()
          else:
            src_series[source_string] = np.nansum(np.nansum(np.ma.array(src_grid, mask=mask_3d),axis=2),axis=1).tolist()
    src_series['totalIN'] = np.nansum(np.nansum(np.ma.array(tot_in, mask=mask_3d),axis=2),axis=1).tolist()
    return src_series

def all_inputs_to_table(params):

    src_series = all_inputs_to_dict(params)

    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'sources_'+get_river_name(params)+'.csv')
    dict_to_csv(filename, src_series)

def all_atm_exch_to_dict(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")

    mask_3d, mask_2d, dummy_name, dummy_nc, dum_asc,  modeldat_startindex, modeldat_endindex, wbody_startindex, wbody_endindex, waterbodyid = \
        make_3d_mask(params, species, folder, 'climate.asc')

    waterbodyid_grid = waterbodyid['waterbodyid'][wbody_startindex:wbody_endindex,:,:]
    waterbodyid.close()

    atm_exch_series = dict()
    atm_exch_series["time"] = manip.convert_numdate2year(dummy_nc['time'][:], dummy_nc['time'].units)[modeldat_startindex:modeldat_endindex]
    for specie in species:
      atm_exch_files = sorted(directory.get_files_with_str(folder, "atmospheric_exchange_"+specie.get_name().upper()+'*.nc'))
      if len(atm_exch_files):
        sub_spec_atm_exch = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
        tot_spec_atm_exch = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
      for ifn in range(len(atm_exch_files)):
        fn = atm_exch_files[ifn]
        nc = Dataset(fn, 'r')
        grid_3d = nc["atmospheric_exchange_"+specie.get_name().upper()][modeldat_startindex:modeldat_endindex,:,:]
        nc.close()
        mask_3d = np.broadcast_to(mask_2d, grid_3d.shape)

        # store individual orders
        if (ifn < len(atm_exch_files)-mainstream_id) or (ifn > len(atm_exch_files)-mainstream_id):
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_order'+str(fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()
        else: # mask the non-outlet lake/reservoirs gridcells in main stream order
          grid_3d_dum = copy.deepcopy(grid_3d)
          grid_3d_lakes = copy.deepcopy(grid_3d)
          grid_3d_reservoirs = copy.deepcopy(grid_3d)
          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
          # store main stream without lakes and reservoirs
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_river'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_order'+str(fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

        # store sum of subgrid orders
        if (ifn < len(atm_exch_files)-mainstream_id):
          sub_spec_atm_exch = np.add(grid_3d, sub_spec_atm_exch)
        elif (ifn == len(atm_exch_files)-mainstream_id):
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_subgrid'] = np.nansum(np.nansum(np.ma.array(sub_spec_atm_exch, mask=mask_3d),axis=2),axis=1).tolist()

        if "atmospheric_exchange_"+specie.get_name().upper()+'_subgrid' in atm_exch_series.keys() and (ifn == len(atm_exch_files)-mainstream_id):
          mask_file = os.path.join(params.water_inputdir, 'channel_width.asc')
          mask_small_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='LE', mask_type='np_grid')
          mask_small_streams_3d_dum = np.broadcast_to(mask_small_streams_2d, grid_3d.shape)
          mask_small_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
          mask_small_streams_3d[:,:,:] = True
          mask_small_streams_3d[np.logical_and(mask_3d==False, mask_small_streams_3d_dum[:,:,:]==False)] = False
          mask_major_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='GT', mask_type='np_grid')
          mask_major_streams_3d_dum = np.broadcast_to(mask_major_streams_2d , grid_3d.shape)
          mask_major_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
          mask_major_streams_3d[:,:,:] = True
          mask_major_streams_3d[np.logical_and(mask_3d==False, mask_major_streams_3d_dum[:,:,:]==False)] = False
          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_smallstreams'] = \
              (np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1)+np.array(atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_subgrid'])).tolist()
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_majorstreams'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist()
          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]

        # store sum of cell
        tot_spec_atm_exch = np.add(grid_3d, tot_spec_atm_exch)

        # store lakes
        if (ifn==len(atm_exch_files)-mainstream_id):
          lakesmask_3d = np.zeros(grid_3d.shape, dtype=bool)
          lakesmask_3d[:,:,:] = True
          lakesmask_3d[np.where(np.logical_and(waterbodyid_grid[:,:,:]>=10000, mask_3d[:,:,:]==False))] = False
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_lakes']  = np.nansum(np.nansum(np.ma.array(grid_3d_lakes, mask=lakesmask_3d),axis=2),axis=1).tolist()

        # store reservoirs
        if (ifn==len(atm_exch_files)-mainstream_id):
          reservoirsmask_3d = np.zeros(grid_3d.shape, dtype=bool)
          reservoirsmask_3d[:,:,:] = True
          reservoirsmask_3d_dum = np.zeros(grid_3d.shape, dtype=bool)
          reservoirsmask_3d_dum[:,:,:] = True

          reservoirsmask_3d_dum[np.where(np.logical_and(waterbodyid_grid[:,:,:]>1, waterbodyid_grid[:,:,:]<10000))] = False
          reservoirsmask_3d[np.where(np.logical_and(reservoirsmask_3d_dum[:,:,:]==False, mask_3d[:,:,:]==False))] = False
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_reservoirs']  = np.nansum(np.nansum(np.ma.array(grid_3d_reservoirs, mask=reservoirsmask_3d),axis=2),axis=1).tolist()

      if len(atm_exch_files)>0:
        atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_total'] = np.nansum(np.nansum(np.ma.array(tot_spec_atm_exch, mask=mask_3d),axis=2),axis=1).tolist()
    return atm_exch_series

def all_atm_exch_to_table(params):
    atm_exch_series = all_atm_exch_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'atm_exch_'+get_river_name(params)+'.csv')
    dict_to_csv(filename, atm_exch_series)

def all_exports_to_dict(params,add_color=False):
    # setattr(params, 'ldebug', True)
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)
    folder = os.path.join(params.outputdir, '..', "BUDGET", "subgrid")

    # mouthmask_3d = np.broadcast_to(mouthmask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()

    mouthmask_3d, mask_2d, dummy_name, dummy_nc, dum_asc,  modeldat_startindex, modeldat_endindex, wbody_startindex, wbody_endindex, waterbodyid = \
        make_3d_mask(params, species, folder, 'rivermouth.asc')

    export_series = dict()
    export_series["time"] = manip.convert_numdate2year(dummy_nc['time'][modeldat_startindex:modeldat_endindex], dummy_nc['time'].units)
    tot_export = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))

    for specie in species:
      all_export = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
      export_files = directory.get_files_with_str(folder, specie.get_name()+"loadOUT_order6*")
      for fn in export_files:
        nc = Dataset(fn, 'r')
        grid_3d = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
        grid_3d[np.where(nc[specie.get_name()+"loadOUT"][modeldat_startindex:modeldat_endindex,:,:]<1e12)] = \
            nc[specie.get_name()+"loadOUT"][modeldat_startindex:modeldat_endindex,:,:][np.where(nc[specie.get_name()+"loadOUT"][modeldat_startindex:modeldat_endindex,:,:]<1e12)]
        nc.close()
        all_export = np.add(grid_3d, all_export)
        if (not specie.get_name().lower()=='alk') and (not 'tss' in specie.get_name().lower()) and (not 'pim' in specie.get_name().lower()):
          tot_export = np.add(grid_3d, tot_export)
      export_series[specie.get_name().upper()+"loadOUT"] = np.nansum(np.nansum(np.ma.array(all_export, mask=mouthmask_3d),axis=2),axis=1).tolist()
      if add_color:
          export_series[specie.get_name().upper()+"loadOUT"] = [np.nansum(np.nansum(np.ma.array(all_export, mask=mouthmask_3d),axis=2),axis=1).tolist(), specie.get_val('color')]
      else:
          export_series[specie.get_name().upper()+"loadOUT"] = np.nansum(np.nansum(np.ma.array(all_export, mask=mouthmask_3d),axis=2),axis=1).tolist()

    export_series["total_loadOUT"] = np.nansum(np.nansum(np.ma.array(tot_export, mask=mouthmask_3d),axis=2),axis=1).tolist()
    return export_series

def all_exports_to_table(params):
    export_series = all_exports_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'exports_'+get_river_name(params)+'.csv')
    dict_to_csv(filename, export_series)


def all_budget_to_dict(params):
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    src_filename = os.path.join(folder, 'sources_'+get_river_name(params)+'.csv')
    atm_filename = os.path.join(folder, 'atm_exch_'+get_river_name(params)+'.csv')
    exp_filename = os.path.join(folder, 'exports_'+get_river_name(params)+'.csv')

    src_series = csv_to_dict(src_filename)
    atm_series = csv_to_dict(atm_filename)
    exp_series = csv_to_dict(exp_filename)


    time = src_series['time'][:]

    dellist = list()
    for key in atm_series.keys():
        if not ('total' in key):
            dellist.append(key)
    for key in dellist:
            del atm_series[key]

    budget_series = dict()
    budget_series['time'] = time
    budget_series['input'] = (np.array(src_series['totalIN'])).tolist()
    budget_series['export'] = (-1*np.array(exp_series['total_loadOUT'])).tolist()
    budget_series['emission'] = (-1*np.array(atm_series["atmospheric_exchange_DIC_total"])).tolist()
    budget_series['retention'] = \
        ((np.array(budget_series['input'])+\
          np.array(budget_series['emission'])+\
          np.array(budget_series['export']))).tolist()
    return budget_series

def all_budget_to_table(params):
    budget_series = all_budget_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'total_budget_'+get_river_name(params)+'.csv')
    dict_to_csv(filename, budget_series)

def all_fluxes_to_dict(params):
    def list_append(lst, item):
      lst.append(item)
      return lst

    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)
    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    # debugprint(params,'print(folder)')
    # debugprint(params,print(folder))

    # mouthmask_3d = np.broadcast_to(mouthmask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()

    mask_3d, mask_2d, dummy_name, dummy_nc, dum_asc,  modeldat_startindex, modeldat_endindex, wbody_startindex, wbody_endindex, waterbodyid = \
        make_3d_mask(params, species, folder, 'climate.asc')

    waterbodyid_grid = waterbodyid['waterbodyid'][wbody_startindex:wbody_endindex,:,:]
    waterbodyid.close()


    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1


    flux_series = dict()
    for specie in species:
      print(specie.name)
      flux_series[specie.get_val('name')] = dict()
      flux_series[specie.get_val('name')]["time"] = manip.convert_numdate2year(dummy_nc['time'][:], dummy_nc['time'].units)[modeldat_startindex:modeldat_endindex]
      stoich = reactions.specie_dy(proc,specie.get_name())
      for iproc in range(len(proc)):
        if stoich[iproc]!=0:
          proc_fns = sorted(directory.get_files_with_str(folder, proc[iproc].get_val("name")+"_order*.nc"))
          sub_spec_flux = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
          tot_spec_flux = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
          for ifn in range(len(proc_fns)):
            proc_fn = proc_fns[ifn]
            proc_nc = Dataset(proc_fn, 'r')
            grid_3d = proc_nc[proc[iproc].get_val("name")][modeldat_startindex:modeldat_endindex,:,:]*stoich[iproc]
            proc_nc.close()

            # store individual orders
            if (ifn < len(proc_fns)-mainstream_id) or (ifn > len(proc_fns)-mainstream_id):
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_order'+str(proc_fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()
            elif (ifn == len(proc_fns)-mainstream_id): # mask the non-outlet lake/reservoirs gridcells in main stream order
              grid_3d_dum = copy.deepcopy(grid_3d)
              grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
			  # store main stream without lakes and reservoirs
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_river'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

              grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_order'+str(proc_fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

            # store sum of subgrid orders
            if (ifn < len(proc_fns)-mainstream_id):
              sub_spec_flux = np.add(grid_3d, sub_spec_flux)
            elif (ifn == len(proc_fns)-mainstream_id):
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_subgrid'] = np.nansum(np.nansum(np.ma.array(sub_spec_flux, mask=mask_3d),axis=2),axis=1).tolist()

            if proc[iproc].get_val("name")+'_subgrid' in flux_series[specie.get_val('name')].keys() and (ifn == len(proc_fns)-mainstream_id):
              mask_file = os.path.join(params.water_inputdir, 'channel_width.asc')
              mask_small_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='LE', mask_type='np_grid')
              mask_small_streams_3d_dum = np.broadcast_to(mask_small_streams_2d, grid_3d.shape)
              mask_small_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
              mask_small_streams_3d[:,:,:] = True
              mask_small_streams_3d[np.logical_and(mask_3d==False, mask_small_streams_3d_dum[:,:,:]==False)] = False
              mask_major_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='GT', mask_type='np_grid')
              mask_major_streams_3d_dum = np.broadcast_to(mask_major_streams_2d , grid_3d.shape)
              mask_major_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
              mask_major_streams_3d[:,:,:] = True
              mask_major_streams_3d[np.logical_and(mask_3d==False, mask_major_streams_3d_dum[:,:,:]==False)] = False
              grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_smallstreams'] = \
              (np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1)+np.array(flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_subgrid'])).tolist()
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_majorstreams'] = \
              np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist()
              grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]

            # store sum of cell
            tot_spec_flux = np.add(grid_3d, tot_spec_flux)

            # store lakes
            if (ifn==len(proc_fns)-mainstream_id):
              lakesmask_3d = np.zeros(grid_3d.shape, dtype=bool)
              lakesmask_3d[:,:,:] = True
              lakesmask_3d[np.where(np.logical_and(waterbodyid_grid[:,:,:]>=10000, mask_3d[:,:,:]==False))] = False
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_lakes']  = np.nansum(np.nansum(np.ma.array(grid_3d, mask=lakesmask_3d),axis=2),axis=1).tolist()

            # store reservoirs
            if (ifn==len(proc_fns)-mainstream_id):
              reservoirsmask_3d = np.zeros(grid_3d.shape, dtype=bool)
              reservoirsmask_3d[:,:,:] = True
              reservoirsmask_3d_dum = np.zeros(grid_3d.shape, dtype=bool)
              reservoirsmask_3d_dum[:,:,:] = True

              reservoirsmask_3d_dum[np.where(np.logical_and(waterbodyid_grid[:,:,:]>1, waterbodyid_grid[:,:,:]<10000))] = False
              reservoirsmask_3d[np.where(np.logical_and(reservoirsmask_3d_dum[:,:,:]==False, mask_3d[:,:,:]==False))] = False
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_reservoirs']  = np.nansum(np.nansum(np.ma.array(grid_3d, mask=reservoirsmask_3d),axis=2),axis=1).tolist()

            # store floodplains
            if 'order7' in proc_fns[ifn]:
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_floodplains/wetlands']  = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()


          flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_totflux'] = np.nansum(np.nansum(np.ma.array(tot_spec_flux, mask=mask_3d),axis=2),axis=1).tolist()

        src_dict = csv_to_dict(os.path.join(params.outputdir, '..', 'ANALYSIS',  "tables", "sources_"+get_river_name(params)+'.csv'))
        if (specie.get_name()+'_srcloadIN') in list(src_dict.keys()):
          if len(src_dict[specie.get_name()+'_srcloadIN']) < len (flux_series[specie.get_val('name')]["time"]):
            flux_series[specie.get_val('name')][specie.get_name()+'_srcflux'] = list_append(src_dict[specie.get_name()+'_srcloadIN'], None)
          else:
            flux_series[specie.get_val('name')][specie.get_name()+'_srcflux'] = src_dict[specie.get_name()+'_srcloadIN']
        else:
          flux_series[specie.get_val('name')][specie.get_name()+'_srcflux'] = [0]*len(flux_series[specie.get_val('name')]["time"])

        exp_fns = directory.get_files_with_str(folder, specie.get_name()+"loadOUT_order6.nc")
        tot_exp = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
        for ifn in range(len(exp_fns)):
            exp_fn = exp_fns[ifn]
            exp_nc = Dataset(exp_fn, 'r')
            grid_3d = exp_nc[specie.get_name()+"loadOUT"][modeldat_startindex:modeldat_endindex,:,:]*-1
            exp_nc.close()
            tot_exp = np.add(grid_3d, tot_exp)
        flux_series[specie.get_val('name')][specie.get_name()+'_expflux'] = np.nansum(np.nansum(np.ma.array(tot_exp, mask=mask_3d),axis=2),axis=1) #mask_3d was mouthmask_3d
      flux_series[specie.get_val('name')]['budget'] = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,0,0]))
      for flux in flux_series[specie.get_val('name')]:
        #debugprint(params,specie.get_val('name'))
        #debugprint(params,flux)
        if 'flux' in flux:
          flux_series[specie.get_val('name')]['budget'] += flux_series[specie.get_val('name')][flux]
    dummy_nc.close()
    return flux_series

def all_fluxes_to_table(params):
    flux_series = all_fluxes_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))
    for speciekey, speciedict in flux_series.items():
      filename = os.path.join(folder, speciekey+'_biogeochemistry_'+get_river_name(params)+'.csv')
      dict_to_csv(filename, speciedict)

def all_sec_to_dict(params):
    folder = os.path.join(params.outputdir, '..', 'STREAM_ENV_CONDITIONS')
    filename_list = directory.get_files_with_str(folder, "*.nc")
    arguments = list()
    for filename in filename_list:
      if (('vol' in filename) or ('area' in filename)) and (not 'dvoldt' in filename):
        arguments.append(os.path.splitext(os.path.basename(filename))[0])

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    define_subgrid_streamorder.define_subgrid_streamorder_constant(params)

    folder = os.path.join(params.outputdir, "..", "STREAM_ENV_CONDITIONS", "subgrid")


    species,sources,proc,params_local = read_parameter.readfile(params.species_ini) # NOT NEEDED, but here 'species' needed for make_3d_mask

    mask_3d, mask_2d, dummy_name, dummy_nc, dum_asc,  modeldat_startindex, modeldat_endindex, wbody_startindex, wbody_endindex, waterbodyid = \
        make_3d_mask(params, species, folder, 'climate.asc')

    waterbodyid_grid = waterbodyid['waterbodyid'][wbody_startindex:wbody_endindex,:,:]
    waterbodyid.close()

    sec_series = dict()
    for arg in arguments:
      sec_series[arg] = dict()
      sec_series[arg]["time"] = manip.convert_numdate2year(dummy_nc['time'][:], dummy_nc['time'].units)[modeldat_startindex:modeldat_endindex]
      arg_fns = sorted(directory.get_files_with_str(folder, arg+"_order*.nc"))
      sub_arg = [0]*(modeldat_endindex-modeldat_startindex)
      tot_arg = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
      for ifn in range(len(arg_fns)):
            arg_fn = arg_fns[ifn]
            arg_nc = Dataset(arg_fn, 'r')
            grid_3d = arg_nc[arg][modeldat_startindex:modeldat_endindex,:,:]
            # store individual orders
            if (ifn < len(arg_fns)-mainstream_id) or (ifn > len(arg_fns)-mainstream_id):
              sec_series[arg][arg+'_order'+str(arg_fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()
            elif (ifn == len(arg_fns)-mainstream_id): # mask the non-outlet lake/reservoirs gridcells in main stream order
              grid_3d_dum = copy.deepcopy(grid_3d)
              grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
			  # store main stream without lakes and reservoirs
              sec_series[arg][arg+'_river'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

              grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]
              sec_series[arg][arg+'_order'+str(arg_fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()


            # store sum of subgrid orders
            if (ifn < len(arg_fns)-mainstream_id):
              sub_arg = np.add(sub_arg, np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1)*params.number_of_rivers[ifn])
            elif (ifn == len(arg_fns)-mainstream_id):
              sec_series[arg][arg+'_subgrid'] = sub_arg

            if (arg+'_subgrid' in sec_series[arg].keys()) and (ifn == len(arg_fns)-mainstream_id):
                mask_file = os.path.join(params.water_inputdir, 'channel_width.asc')
                mask_small_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='LE', mask_type='np_grid')
                mask_small_streams_3d_dum = np.broadcast_to(mask_small_streams_2d, grid_3d.shape)
                mask_small_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
                mask_small_streams_3d[:,:,:] = True
                mask_small_streams_3d[np.logical_and(mask_3d==False, mask_small_streams_3d_dum[:,:,:]==False)] = False
                mask_major_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='GT', mask_type='np_grid')
                mask_major_streams_3d_dum = np.broadcast_to(mask_major_streams_2d , grid_3d.shape)
                mask_major_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
                mask_major_streams_3d[:,:,:] = True
                mask_major_streams_3d[np.logical_and(mask_3d==False, mask_major_streams_3d_dum[:,:,:]==False)] = False
                grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
                sec_series[arg][arg+'_smallstreams'] = (np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1)+np.array(sec_series[arg][arg+'_subgrid'])).tolist()
                sec_series[arg][arg+'_majorstreams'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist()
                grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]

            # store sum of cell
            tot_arg = np.add(grid_3d, tot_arg)

            # store lakes
            if (ifn==len(arg_fns)-mainstream_id):
              lakesmask_3d = np.zeros(grid_3d.shape, dtype=bool)
              lakesmask_3d[:,:,:] = True
              lakesmask_3d[np.where(np.logical_and(waterbodyid_grid[:,:,:]>=10000, mask_3d[:,:,:]==False))] = False
              sec_series[arg][arg+'_lakes']  = np.nansum(np.nansum(np.ma.array(grid_3d, mask=lakesmask_3d),axis=2),axis=1).tolist()

            # store reservoirs
            if (ifn==len(arg_fns)-mainstream_id):
              reservoirsmask_3d = np.zeros(grid_3d.shape, dtype=bool)
              reservoirsmask_3d[:,:,:] = True
              reservoirsmask_3d_dum = np.zeros(grid_3d.shape, dtype=bool)
              reservoirsmask_3d_dum[:,:,:] = True

              reservoirsmask_3d_dum[np.where(np.logical_and(waterbodyid_grid[:,:,:]>1, waterbodyid_grid[:,:,:]<10000))] = False
              reservoirsmask_3d[np.where(np.logical_and(reservoirsmask_3d_dum[:,:,:]==False, mask_3d[:,:,:]==False))] = False
              sec_series[arg][arg+'_reservoirs']  = np.nansum(np.nansum(np.ma.array(grid_3d, mask=reservoirsmask_3d),axis=2),axis=1).tolist()

            # store floodplains
            if 'order7' in arg_fns[ifn]:
              sec_series[arg][arg+'_floodplains/wetlands']  = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()
              sec_series[arg].pop(arg+'_order7', None)


      sec_series[arg][arg+'_tot'] = np.nansum(np.nansum(np.ma.array(tot_arg, mask=mask_3d),axis=2),axis=1).tolist()
    dummy_nc.close()
    return sec_series

def all_stream_env_conditions_to_table(params):
    sec_series = all_sec_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))
    for key, dictionary in sec_series.items():
      filename = os.path.join(folder, key+'_sec_'+get_river_name(params)+'.csv')
      dict_to_csv(filename, dictionary)

# def map_net_budget(params):
#     species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
#     make_index.species(params,species,proc)
#     basinid = pickle.load(open(os.path.join(params.outputdir,'mask_asc.pkl'),'rb'))
#     folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
#     dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
#     mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
#     climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
#     climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
#     mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
#     mask_2d[:,:] = True
#     mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False
#     try:
#         all_budget
#     except NameError:
#         proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

#         dummy_nc = Dataset(proclist[0], 'r')
#         dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]

#         specie_budget = np.zeros(np.shape(dummy_nc[dummy_name]))
#         all_budget = np.zeros(np.shape(dummy_nc[dummy_name]))
#     else:
#         pass

#     for specie in species:
#         stoich = reactions.specie_dy(proc,specie.get_name())
#         specie_budget = np.zeros(np.shape(dummy_nc[dummy_name]))
#         for iproc in range(len(proc)):
#             proc_fn = os.path.join(folder, proc[iproc].get_val("name")+"_order6.nc")
#             proc_nc = Dataset(proc_fn, 'r')
#             grid_3d = proc_nc[proc[iproc].get_val("name")][:,:,:]*stoich[iproc]
#             specie_budget = np.add(specie_budget, grid_3d)
#         # also account for advection terms
#         adv_terms = ['upstrloadIN','loadOUT','srcloadIN']
#         adv_pos = [1,-1, 1]
#         for i in range(len(adv_terms)):
#             proc_fn = os.path.join(folder, specie.get_name()+'_'+adv_terms[i]+"_order6.nc")
#             proc_nc = Dataset(proc_fn, 'r')
#             grid_3d = proc_nc[specie.get_name()+'_'+adv_terms[i]][:,:,:]*adv_pos[i]
#             specie_budget = np.add(specie_budget, grid_3d)
#         mask_3d = np.broadcast_to(mask_2d, specie_budget.shape)
#         specie_budget = np.ma.array(specie_budget, mask=mask_3d)
#         if ((specie.get_name().lower()!='alk') and (not 'tss' in specie.get_name().lower())):
#           all_budget += np.add(all_budget, specie_budget)


#         specie_budget_nc = Dataset(os.path.join(folder, specie.get_name()+"_budget.nc"), 'w')
#         output_conversion.init_ncdata(folder, specie_budget_nc, specie.get_name()+"_budget", basinid, unit='Mmol/yr', long_name=specie.get_name()+"_budget")
#         for itime in range(len(dummy_nc['time'])):
#             manip.add_grid_time(specie_budget_nc, specie.get_name()+"_budget", specie_budget[itime,:,:], dummy_nc['time'][itime])

#     all_budget = np.ma.array(all_budget, mask=mask_3d)
#     all_budget_nc = Dataset(os.path.join(folder, "all_budget.nc"), 'w')
#     output_conversion.init_ncdata(folder, all_budget_nc, "all_budget", basinid, unit='Mmol/yr', long_name="all_budget")
#     for itime in range(len(dummy_nc['time'])):
#         manip.add_grid_time(all_budget_nc, "all_budget", all_budget[itime,:,:], dummy_nc['time'][itime])

if __name__ == "__main__":
    # Set the general path for the own python modules

    # import general_path
    # Parse command-line arguments and set parameters for script
    #try:
    # Parse command-line arguments and set parameters for script
    # Startup logging and runtime start
    params,log,s = general_startup.general_startup(sys.argv)
    #except SystemExit:
    #  raise MyError("Error has occured in the reading of the commandline options.")

    do(params)