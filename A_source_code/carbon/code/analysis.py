from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import imageio
import csv
import copy
import pandas as pd

import os
import sys

import general_path 
    
import read_parameter

import ascraster
import cmd_options_dgnm
import directory
import get_dynamic_gridinfo
import general_startup
import interpolate_list
import make_index_species
import make_outlakes_dict
#import make_plots
import make_mask
import output_conversion
import pickle
import pointer_class
import reactions
import specie
import define_subgrid_streamorder

import manip

global_colors = cm.Set1(np.linspace(0,1,8))
table_start_year=1950
table_end_year=2000

def do(params):

    ## make all tables

    all_inputs_to_table(params)
    #all_atm_exch_to_table(params)
    #all_atm_exch_to_climate_tables(params)
    #all_burial_to_table(params)
    #all_burial_to_climate_tables(params)
    #all_exports_to_table(params)
    #all_budget_to_table(params)
    #all_fluxes_to_table(params)

    ## make latitudinal plots (last timestep)
    #if (params.maskid==0) and (params.mask_bool_operator=='GT'):
      #all_sec_per_latitude(params)
      #all_inputs_per_latitude(params)
      #all_atm_exch_per_latitude(params)
      #all_fluxes_per_latitude(params)
   
    #all_fluxes_to_multiplots(params)
    #all_speciefluxes_to_multiplots(params)

    #all_fluxes_to_stackserieplots(params)
    all_fluxes_to_stackserieplots_TgPerYear(params)
    #all_fluxes_to_maps(params)
    all_budget_to_stackserieplots_Tg_year(params)
    #map_net_budget(params)

    all_stream_env_conditions_to_table(params)
    if params.maskid==0 and params.mask_bool_operator=='GT':
      all_stream_env_conditions_to_climate_tables(params)
      all_fluxes_to_climate_tables(params)
    conv_all_tables_to_Tg(params)


    ## make flux tables per LOAC component
    #all_fluxes_to_smallstreams_to_table(params)
    #all_fluxes_to_rivers_to_table(params)
    #all_fluxes_to_lakes_to_table(params)
    #all_fluxes_to_reservoirs_to_table(params)
    #all_fluxes_to_floodplains_to_table(params)

    ## make all pieplots
    #all_inputs_to_pieplots(params)
    #all_exports_to_pieplots(params)
    #all_fluxes_to_pieplots(params)
    #ll_fluxes_per_wbtype_to_pieplots(params)
    #all_fluxes_per_basin_to_pieplots(params)
    #all_fluxes_per_climate_to_pieplots(params)

    #all_inputs_to_stackserieplots(params)
    #all_atm_exch_to_plots(params)
    #all_atm_exch_to_stackserieplots_MmolPerTimestep(params)
    #all_atm_exch_to_stackserieplots_MmolPerYear(params)
    #all_atm_exch_to_stackserieplots(params)

    #all_exports_to_stackserieplots(params)
	    
def get_river_name(params):
  if (not 'country' in params.file_mask) and (params.mask_bool_operator=='EQ'):
    if params.maskid==99:
      rivername = "Rhine"
    elif params.maskid==1:
      rivername = 'Amazon'
    elif params.maskid==2:
      rivername = 'Congo'
    elif params.maskid==3:
      rivername = 'Caspian Sea'
    elif params.maskid==4:
      rivername = 'Mississippi'
    elif params.maskid==5:
      rivername = 'Nile'
    elif params.maskid==10:
      rivername = 'Yangtze'
    elif params.maskid==11:
      rivername = 'Amur'
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
'''
def make_colors(params):
    color_dict = dict()
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)

    i = 0
    for spec in species:
      color_dict[key] = global_colors[i]
      i+=1
    return color_dict
'''

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

def make_time_indices(params):
    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1  
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1       
    return modeldat_startindex, modeldat_endindex, all_dat_startindex, all_dat_endindex

def make_3d_mask(mask_2d, modeldat_startindex, modeldat_endindex, dummy_nc, dummy_name):
    return np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape) 

def dict_to_csv(filename, mydict):
    pd.DataFrame(mydict).T.reset_index().to_csv(filename, header=False, index=False)

def csv_to_dict(filename):
    dummy_dict = pd.read_csv(filename, index_col=0, header=None).to_dict("split")
    mydict = dict(zip(dummy_dict['index'], dummy_dict['data']))
    return mydict

def conv_all_tables_to_Tg(params):
    table_dir = os.path.join(params.outputdir, "..", "ANALYSIS", "tables")
    table_list = directory.get_files_with_str(table_dir, "burial*.csv", exclude=['validation', 'Tg'])
    table_list.extend(directory.get_files_with_str(table_dir, "atm_exch*.csv", exclude=['validation', 'Tg']))
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

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    #mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     
    mask_2d[np.where(mask_2d_dum[:,:]==False)] = False     

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    src_folder = params.load_inputdir
    
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")
    
    dummy_nc = Dataset(proclist[-1], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[-1]))[0][:-7]

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)   

    #print("all_dat_startindex ",all_dat_startindex)
    #print(" ",)
    #print(" ",)
    #print(" ",)

    src_series = dict()
    src_series["time"] = manip.convert_numdate2year(dummy_nc['time'][modeldat_startindex:modeldat_endindex], dummy_nc['time'].units) 

    tot_in = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))

    for specie in species:
      src_files = sorted(directory.get_files_with_str(folder, "atmospheric_exchange_"+specie.get_name().upper()+'*.nc'))    

    ## this is still a bit too specific for carbon; to be formulated more generic
    for source in sources:
      src_nc = Dataset(os.path.join(src_folder, source.get_val('name')+'.nc'), 'r')
      for attrib in source.get_attrib():
          if 'fr_' in attrib:
            if (not 'TSS' in attrib) and (not 'PIM' in attrib):
              fraction = getattr(source, attrib)
            elif ('tss' in source.get_val('name').lower() and ('TSS' in attrib)) or ('pim' in source.get_val('name').lower() and ('PIM' in attrib)):
              fraction = 1
      src_grid = src_nc[source.get_val('name')][all_dat_startindex:all_dat_endindex,:,:]*params.outputtime*fraction
      src_series[source.get_val('name')] = np.nansum(np.nansum(np.ma.array(src_grid, mask=mask_3d),axis=2),axis=1).tolist()
      for attrib in source.get_attrib():
        if 'fr_' in attrib:
          if (not 'TSS' in attrib) and (not 'ALK' in attrib) and (not 'PIM' in attrib):
            tot_in = np.add(src_grid, tot_in)
      for attrib in source.get_attrib():
        if 'fr_' in attrib:
          src_grid = src_nc[source.get_val('name')][all_dat_startindex:all_dat_endindex,:,:]*params.outputtime*fraction
          
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

def all_inputs_to_pieplots(params):
    src_series_dum = all_inputs_to_dict(params,add_color=True)
    src_series = copy.deepcopy(src_series_dum)
    directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))
    for key in src_series_dum:
      if not 'srcloadIN' in key:
        del src_series[key]
 
    src_values = np.array([])
    src_labels = np.array([])
    for key,value_list in src_series.items():
        value = np.array(value_list)
        if (not key == 'time') and (not 'total' in key) and (not 'tss' in key.lower()) and not (not 'pim' in key.lower()) and (value.size>0) and (not 'alk' in key.lower()) and (not 'benth' in key.lower()):
          src_values = np.append(src_values, np.mean(value[int(value.size*0.5):]))
          src_labels = np.append(src_labels, key+' ['+'{0:.2f}'.format((np.mean(value[int(value.size*0.5):])/np.mean(src_series_dum['totalIN'][int(value.size*0.5):]))*100)+"%]")
    annotations = ['Results from: \n'+params.outputdir]
    outname = os.path.join(params.outputdir, '..', 'ANALYSIS', "plots", get_river_name(params)+'_PIECHART_sources.png')
    make_plots.piechart(src_values, src_labels, outfilename=outname, title='Carbon inputs to the '+get_river_name(params), annotations=annotations, colors='Blues')

def all_inputs_to_stackserieplots(params):      
    filename = os.path.join(params.outputdir, '..', 'ANALYSIS', "tables", 'sources_'+get_river_name(params)+'.csv')
    src_series = csv_to_dict(filename)
    directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))
    time = src_series['time']
    dellist = list()
    conv = 12*1e-6*(1/params.outputtime)
    for key in src_series.keys():
        if ('tot' in key) or ('lakes' in key) or ('reservoir' in key) or ('time' in key) or ('TSS' in key) or ('ALK' in key) or ('PIM' in key) or not ('srcloadIN' in key):
            dellist.append(key)
        else:
          src_series[key]=(np.array(src_series[key])*conv).tolist()
    for key in dellist:
            del src_series[key]

    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'total_budget_'+get_river_name(params)+'.csv')
    budget_series = csv_to_dict(filename)
    y_min, y_max = 0, 1.1*max(budget_series['input'])*conv
    color_dict = make_colors(params)

    outname = os.path.join(params.outputdir, '..', 'ANALYSIS', "plots", get_river_name(params)+'_TIMESERIES_sources.png')	
    make_plots.stacked_linechart(time, src_series, \
                                           xlab='time [years]', ylab='flux [Tg/yr]', \
                                           outfilename=outname, \
                                           title='Delivery of carbon forms to '+get_river_name(params), 
                                           color='mixed', color_dict=color_dict, ylim=(0, y_max), xlim=(budget_series['time'][0]-3,budget_series['time'][-1]+3))

def all_inputs_per_latitude(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    src_folder = params.load_inputdir
    
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")
    
    dummy_nc = Dataset(proclist[-1], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[-1]))[0][:-7]

    conv = 12*1e-6*(1/params.outputtime)

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)   

    for specie in species:
      src_files = sorted(directory.get_files_with_str(folder, "atmospheric_exchange_"+specie.get_name().upper()+'*.nc'))    

    lat_array = np.linspace(89.75,-89.75, 360)
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "plots"))
    

    ## this is still a bit too specific for carbon; to be formulated more generic
    for source in sources:
      src_nc = Dataset(os.path.join(src_folder, source.get_val('name')+'.nc'), 'r')
      for attrib in source.get_attrib():
          if 'fr_' in attrib:
            if (not 'TSS' in attrib) and (not 'PIM' in attrib):
              fraction = getattr(source, attrib)
            elif ('tss' in source.get_val('name').lower() and ('TSS' in attrib)) or ('pim' in source.get_val('name').lower() and ('PIM' in attrib)):
              fraction = 1

      for attrib in source.get_attrib():
        if 'fr_' in attrib:
          src_grid = src_nc[source.get_val('name')][-1,:,:]*params.outputtime*fraction*conv
          src_lat_array = np.nansum(np.ma.array(src_grid, mask=mask_2d),axis=1).tolist()
          filename = os.path.join(folder, 'lat_'+source.get_val('name')+'.png')
          make_plots.single_2d(filename, src_lat_array, lat_array, title="None", ylabel="None", xlabel="None")

def all_burial_to_dict(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    # read river basin map
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]   

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  


    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]

    waterbodyid.close()
    endo_waterbodyid.close()

    burial_series = dict()
    burial_series["time"] = manip.convert_numdate2year(dummy_nc['time'][:], dummy_nc['time'].units)[modeldat_startindex:modeldat_endindex]
    waterbodyoutlet.close()
    for specie in species:
      if not 'PIM' in specie.get_name().upper():
        burial_files = sorted(directory.get_files_with_str(folder, "burial_"+specie.get_name()+'*.nc'))
      else:
        burial_files = []
      if len(burial_files):
        sub_spec_burial = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
        tot_spec_burial = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
      for ifn in range(len(burial_files)):
        fn = burial_files[ifn]
        nc = Dataset(fn, 'r')
        grid_3d = nc["burial_"+specie.get_name()][modeldat_startindex:modeldat_endindex,:,:]
        nc.close()
        mask_3d = np.broadcast_to(mask_2d, grid_3d.shape)

        # store individual orders
        if (ifn < len(burial_files)-mainstream_id) or (ifn > len(burial_files)-mainstream_id):
          burial_series["burial_"+specie.get_name().upper()+'_order'+str(fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()
        else: # mask the non-outlet lake/reservoirs gridcells in main stream order
          grid_3d_dum = copy.deepcopy(grid_3d)  
          grid_3d_lakes = copy.deepcopy(grid_3d) 
          grid_3d_reservoirs = copy.deepcopy(grid_3d) 
          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
          # store main stream without lakes and reservoirs
          burial_series["burial_"+specie.get_name().upper()+'_river'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]
          burial_series["burial_"+specie.get_name().upper()+'_order'+str(fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

        # store sum of subgrid orders
        if (ifn < len(burial_files)-mainstream_id):
          sub_spec_burial = np.add(grid_3d, sub_spec_burial)
        elif (ifn == len(burial_files)-mainstream_id):
          burial_series["burial_"+specie.get_name().upper()+'_subgrid'] = np.nansum(np.nansum(np.ma.array(sub_spec_burial, mask=mask_3d),axis=2),axis=1).tolist()

        if "burial_"+specie.get_name().upper()+'_subgrid' in burial_series.keys() and (ifn == len(burial_files)-mainstream_id):
          mask_file = os.path.join(params.water_inputdir, 'channel_width.asc')
          mask_small_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='LE', mask_type='np_grid')
          mask_small_streams_3d_dum = np.broadcast_to(mask_small_streams_2d, grid_3d.shape)
          mask_small_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
          mask_small_streams_3d[:,:,:] = True
          mask_small_streams_3d[np.where(np.logical_and(mask_3d==False, mask_small_streams_3d_dum[:,:,:]==False))] = False
          mask_major_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='GT', mask_type='np_grid')
          mask_major_streams_3d_dum = np.broadcast_to(mask_major_streams_2d , grid_3d.shape)
          mask_major_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
          mask_major_streams_3d[:,:,:] = True
          mask_major_streams_3d[np.where(np.logical_and(mask_3d==False, mask_major_streams_3d_dum[:,:,:]==False))] = False
          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
          burial_series["burial_"+specie.get_name().upper()+'_smallstreams'] = (np.add(np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1),np.array(burial_series["burial_"+specie.get_name().upper()+'_subgrid']))).tolist()
          burial_series["burial_"+specie.get_name().upper()+'_majorstreams'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist() 
          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]

        # store sum of cell
        tot_spec_burial = np.add(grid_3d, tot_spec_burial)

        # store lakes
        if (ifn==len(burial_files)-mainstream_id):
          lakesmask_3d = np.zeros(grid_3d.shape, dtype=bool)
          lakesmask_3d[:,:,:] = True
          lakesmask_3d[np.where(np.logical_and(waterbodyid_grid[:,:,:]>=10000, mask_3d[:,:,:]==False))] = False     
          burial_series["burial_"+specie.get_name().upper()+'_lakes']  = np.nansum(np.nansum(np.ma.array(grid_3d_lakes, mask=lakesmask_3d),axis=2),axis=1).tolist()

        # store reservoirs
        if (ifn==len(burial_files)-mainstream_id):
          reservoirsmask_3d = np.zeros(grid_3d.shape, dtype=bool)
          reservoirsmask_3d[:,:,:] = True
          reservoirsmask_3d_dum = np.zeros(grid_3d.shape, dtype=bool)
          reservoirsmask_3d_dum[:,:,:] = True

          reservoirsmask_3d_dum[np.where(np.logical_and(waterbodyid_grid[:,:,:]>1, waterbodyid_grid[:,:,:]<10000))] = False   
          reservoirsmask_3d[np.where(np.logical_and(reservoirsmask_3d_dum[:,:,:]==False, mask_3d[:,:,:]==False))] = False    
          burial_series["burial_"+specie.get_name().upper()+'_reservoirs']  = np.nansum(np.nansum(np.ma.array(grid_3d_reservoirs, mask=reservoirsmask_3d),axis=2),axis=1).tolist()

    if len(burial_files)>0:
        burial_series["burial_"+specie.get_name().upper()+'_total'] = np.nansum(np.nansum(np.ma.array(tot_spec_burial, mask=mask_3d),axis=2),axis=1).tolist()
    return burial_series

def burial_to_climate_dict(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    # read river basin map
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]   

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d_dum = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d_dum = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  


    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]

    waterbodyid.close()
    endo_waterbodyid.close()
    waterbodyoutlet.close()

    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    tropics_mask_2d = make_mask.do(climatemask_fn, 1, dum_asc, mask_type='np_grid', logical='EQ')
    tropics_mask_3d = np.broadcast_to(tropics_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    arid_mask_2d = make_mask.do(climatemask_fn, 2, dum_asc, mask_type='np_grid', logical='EQ')
    arid_mask_3d = np.broadcast_to(arid_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    temperate_mask_2d = make_mask.do(climatemask_fn, 3, dum_asc, mask_type='np_grid', logical='EQ')
    temperate_mask_3d = np.broadcast_to(temperate_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    continental_mask_2d = make_mask.do(climatemask_fn, 4, dum_asc, mask_type='np_grid', logical='EQ')
    continental_mask_3d = np.broadcast_to(continental_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    polar_mask_2d = make_mask.do(climatemask_fn, 5, dum_asc, mask_type='np_grid', logical='EQ')
    polar_mask_3d = np.broadcast_to(polar_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
 
    name_list = ['tropical', 'arid', 'temperate', 'continental', 'polar']
    climate_mask_list = [tropics_mask_3d, arid_mask_3d, temperate_mask_3d, continental_mask_3d, polar_mask_3d]
    burial_series_list = []
    for climate_name in name_list:
      burial_series = dict()
      burial_series["time"] = manip.convert_numdate2year(dummy_nc['time'][:], dummy_nc['time'].units)[modeldat_startindex:modeldat_endindex]
      mask_3d = np.zeros(mask_3d_dum.shape, dtype=bool)
      mask_3d[:,:,:] = True
      mask_3d[np.logical_and(climate_mask_list[name_list.index(climate_name)][:,:,:]==False, mask_3d_dum[:,:,:]==False)] = False
      for specie in species:
        if not 'PIM' in specie.get_name().upper():
          burial_files = sorted(directory.get_files_with_str(folder, "burial_"+specie.get_name()+'*.nc'))
        else:
          burial_files = []
        if len(burial_files):
          sub_spec_burial = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
          tot_spec_burial = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
        for ifn in range(len(burial_files)):
          fn = burial_files[ifn]
          nc = Dataset(fn, 'r')
          grid_3d = nc["burial_"+specie.get_name()][modeldat_startindex:modeldat_endindex,:,:]
          nc.close()

          # store individual orders
          if (ifn < len(burial_files)-mainstream_id) or (ifn > len(burial_files)-mainstream_id):
            burial_series["burial_"+specie.get_name().upper()+'_order'+str(fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()
          else: # mask the non-outlet lake/reservoirs gridcells in main stream order
            grid_3d_dum = copy.deepcopy(grid_3d)  
            grid_3d_lakes = copy.deepcopy(grid_3d) 
            grid_3d_reservoirs = copy.deepcopy(grid_3d) 
            grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
            # store main stream without lakes and reservoirs
            burial_series["burial_"+specie.get_name().upper()+'_river'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

            grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]
            burial_series["burial_"+specie.get_name().upper()+'_order'+str(fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

          # store sum of subgrid orders
          if (ifn < len(burial_files)-mainstream_id):
            sub_spec_burial = np.add(grid_3d, sub_spec_burial)
          elif (ifn == len(burial_files)-mainstream_id):
            burial_series["burial_"+specie.get_name().upper()+'_subgrid'] = np.nansum(np.nansum(np.ma.array(sub_spec_burial, mask=mask_3d),axis=2),axis=1).tolist()

          if "burial_"+specie.get_name().upper()+'_subgrid' in burial_series.keys() and (ifn == len(burial_files)-mainstream_id):
            mask_file = os.path.join(params.water_inputdir, 'channel_width.asc')
            mask_small_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='LE', mask_type='np_grid')
            mask_small_streams_3d_dum = np.broadcast_to(mask_small_streams_2d, grid_3d.shape)
            mask_small_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
            mask_small_streams_3d[:,:,:] = True
            mask_small_streams_3d[np.where(np.logical_and(mask_3d==False, mask_small_streams_3d_dum[:,:,:]==False))] = False
            mask_major_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='GT', mask_type='np_grid')
            mask_major_streams_3d_dum = np.broadcast_to(mask_major_streams_2d , grid_3d.shape)
            mask_major_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
            mask_major_streams_3d[:,:,:] = True
            mask_major_streams_3d[np.where(np.logical_and(mask_3d==False, mask_major_streams_3d_dum[:,:,:]==False))] = False
            grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
            burial_series["burial_"+specie.get_name().upper()+'_smallstreams'] = (np.add(np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1),np.array(burial_series["burial_"+specie.get_name().upper()+'_subgrid']))).tolist()
            burial_series["burial_"+specie.get_name().upper()+'_majorstreams'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist() 
            grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]

          # store sum of cell
          tot_spec_burial = np.add(grid_3d, tot_spec_burial)

          # store lakes
          if (ifn==len(burial_files)-mainstream_id):
            lakesmask_3d = np.zeros(grid_3d.shape, dtype=bool)
            lakesmask_3d[:,:,:] = True
            lakesmask_3d[np.where(np.logical_and(waterbodyid_grid[:,:,:]>=10000, mask_3d[:,:,:]==False))] = False     
            burial_series["burial_"+specie.get_name().upper()+'_lakes']  = np.nansum(np.nansum(np.ma.array(grid_3d_lakes, mask=lakesmask_3d),axis=2),axis=1).tolist()

          # store reservoirs
          if (ifn==len(burial_files)-mainstream_id):
            reservoirsmask_3d = np.zeros(grid_3d.shape, dtype=bool)
            reservoirsmask_3d[:,:,:] = True
            reservoirsmask_3d_dum = np.zeros(grid_3d.shape, dtype=bool)
            reservoirsmask_3d_dum[:,:,:] = True

            reservoirsmask_3d_dum[np.where(np.logical_and(waterbodyid_grid[:,:,:]>1, waterbodyid_grid[:,:,:]<10000))] = False   
            reservoirsmask_3d[np.where(np.logical_and(reservoirsmask_3d_dum[:,:,:]==False, mask_3d[:,:,:]==False))] = False    
            burial_series["burial_"+specie.get_name().upper()+'_reservoirs']  = np.nansum(np.nansum(np.ma.array(grid_3d_reservoirs, mask=reservoirsmask_3d),axis=2),axis=1).tolist()

        if len(burial_files)>0:
          burial_series["burial_"+specie.get_name().upper()+'_total'] = np.nansum(np.nansum(np.ma.array(tot_spec_burial, mask=mask_3d),axis=2),axis=1).tolist()
      burial_series_list.append(burial_series)
    return burial_series_list

def all_burial_to_table(params):
    burial_series = all_burial_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'burial_'+get_river_name(params)+'.csv')
    dict_to_csv(filename, burial_series)

def all_burial_to_climate_tables(params):
    name_list = ['tropical', 'arid', 'temperate', 'continental', 'polar']
    burial_series_list = burial_to_climate_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))
    i=0
    for burial_series in burial_series_list:
      filename = os.path.join(folder, 'burial_'+name_list[i]+get_river_name(params)+'.csv')
      dict_to_csv(filename, burial_series)
      i+=1


def all_atm_exch_to_dict(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    # read river basin map
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]   

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  


    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]

    waterbodyid.close()
    endo_waterbodyid.close()

    atm_exch_series = dict()
    atm_exch_series["time"] = manip.convert_numdate2year(dummy_nc['time'][:], dummy_nc['time'].units)[modeldat_startindex:modeldat_endindex]
    waterbodyoutlet.close()
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
          mask_small_streams_3d[np.where(np.logical_and(mask_3d==False, mask_small_streams_3d_dum[:,:,:]==False))] = False
          mask_major_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='GT', mask_type='np_grid')
          mask_major_streams_3d_dum = np.broadcast_to(mask_major_streams_2d , grid_3d.shape)
          mask_major_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
          mask_major_streams_3d[:,:,:] = True
          mask_major_streams_3d[np.where(np.logical_and(mask_3d==False, mask_major_streams_3d_dum[:,:,:]==False))] = False
          grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
          atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_smallstreams'] = (np.add(np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1),np.array(atm_exch_series["atmospheric_exchange_"+specie.get_name().upper()+'_subgrid']))).tolist()
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

def atm_exch_to_climate_dict(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    # read river basin map
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]   

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d_dum = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d_dum = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  


    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]

    waterbodyid.close()
    endo_waterbodyid.close()
    waterbodyoutlet.close()

    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    tropics_mask_2d = make_mask.do(climatemask_fn, 1, dum_asc, mask_type='np_grid', logical='EQ')
    tropics_mask_3d = np.broadcast_to(tropics_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    arid_mask_2d = make_mask.do(climatemask_fn, 2, dum_asc, mask_type='np_grid', logical='EQ')
    arid_mask_3d = np.broadcast_to(arid_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    temperate_mask_2d = make_mask.do(climatemask_fn, 3, dum_asc, mask_type='np_grid', logical='EQ')
    temperate_mask_3d = np.broadcast_to(temperate_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    continental_mask_2d = make_mask.do(climatemask_fn, 4, dum_asc, mask_type='np_grid', logical='EQ')
    continental_mask_3d = np.broadcast_to(continental_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    polar_mask_2d = make_mask.do(climatemask_fn, 5, dum_asc, mask_type='np_grid', logical='EQ')
    polar_mask_3d = np.broadcast_to(polar_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
 
    name_list = ['tropical', 'arid', 'temperate', 'continental', 'polar']
    climate_mask_list = [tropics_mask_3d, arid_mask_3d, temperate_mask_3d, continental_mask_3d, polar_mask_3d]
    atm_exch_series_list = []
    for climate_name in name_list:
      atm_exch_series = dict()
      atm_exch_series["time"] = manip.convert_numdate2year(dummy_nc['time'][:], dummy_nc['time'].units)[modeldat_startindex:modeldat_endindex]
      mask_3d = np.zeros(mask_3d_dum.shape, dtype=bool)
      mask_3d[:,:,:] = True
      mask_3d[np.logical_and(climate_mask_list[name_list.index(climate_name)][:,:,:]==False, mask_3d_dum[:,:,:]==False)] = False
      for specie in species:
        atm_exch_files = sorted(directory.get_files_with_str(folder, "atmospheric_exchange_"+specie.get_name()+'*.nc'))
        if len(atm_exch_files):
          sub_spec_atm_exch = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
          tot_spec_atm_exch = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
        for ifn in range(len(atm_exch_files)):
          fn = atm_exch_files[ifn]
          nc = Dataset(fn, 'r')
          grid_3d = nc["atmospheric_exchange_"+specie.get_name()][modeldat_startindex:modeldat_endindex,:,:]
          nc.close()

          # store individual orders
          if (ifn < len(atm_exch_files)-mainstream_id) or (ifn > len(atm_exch_files)-mainstream_id):
            atm_exch_series["atm_exch_"+specie.get_name().upper()+'_order'+str(fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()
          else: # mask the non-outlet lake/reservoirs gridcells in main stream order
            grid_3d_dum = copy.deepcopy(grid_3d)  
            grid_3d_lakes = copy.deepcopy(grid_3d) 
            grid_3d_reservoirs = copy.deepcopy(grid_3d) 
            grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
            # store main stream without lakes and reservoirs
            atm_exch_series["atm_exch_"+specie.get_name().upper()+'_river'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

            grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]
            atm_exch_series["atm_exch_"+specie.get_name().upper()+'_order'+str(fn[-4])] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2),axis=1).tolist()

          # store sum of subgrid orders
          if (ifn < len(atm_exch_files)-mainstream_id):
            sub_spec_atm_exch = np.add(grid_3d, sub_spec_atm_exch)
          elif (ifn == len(atm_exch_files)-mainstream_id):
            atm_exch_series["atm_exch_"+specie.get_name().upper()+'_subgrid'] = np.nansum(np.nansum(np.ma.array(sub_spec_atm_exch, mask=mask_3d),axis=2),axis=1).tolist()

          if "atm_exch_"+specie.get_name().upper()+'_subgrid' in atm_exch_series.keys() and (ifn == len(atm_exch_files)-mainstream_id):
            mask_file = os.path.join(params.water_inputdir, 'channel_width.asc')
            mask_small_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='LE', mask_type='np_grid')
            mask_small_streams_3d_dum = np.broadcast_to(mask_small_streams_2d, grid_3d.shape)
            mask_small_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
            mask_small_streams_3d[:,:,:] = True
            mask_small_streams_3d[np.where(np.logical_and(mask_3d==False, mask_small_streams_3d_dum[:,:,:]==False))] = False
            mask_major_streams_2d = make_mask.do(mask_file, 60, dum_asc, logical='GT', mask_type='np_grid')
            mask_major_streams_3d_dum = np.broadcast_to(mask_major_streams_2d , grid_3d.shape)
            mask_major_streams_3d = np.zeros(mask_3d.shape, dtype=bool)
            mask_major_streams_3d[:,:,:] = True
            mask_major_streams_3d[np.where(np.logical_and(mask_3d==False, mask_major_streams_3d_dum[:,:,:]==False))] = False
            grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
            atm_exch_series["atm_exch_"+specie.get_name().upper()+'_smallstreams'] = (np.add(np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1),np.array(atm_exch_series["atm_exch_"+specie.get_name().upper()+'_subgrid']))).tolist()
            atm_exch_series["atm_exch_"+specie.get_name().upper()+'_majorstreams'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist() 
            grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]

          # store sum of cell
          tot_spec_atm_exch = np.add(grid_3d, tot_spec_atm_exch)

          # store lakes
          if (ifn==len(atm_exch_files)-mainstream_id):
            lakesmask_3d = np.zeros(grid_3d.shape, dtype=bool)
            lakesmask_3d[:,:,:] = True
            lakesmask_3d[np.where(np.logical_and(waterbodyid_grid[:,:,:]>=10000, mask_3d[:,:,:]==False))] = False     
            atm_exch_series["atm_exch_"+specie.get_name().upper()+'_lakes']  = np.nansum(np.nansum(np.ma.array(grid_3d_lakes, mask=lakesmask_3d),axis=2),axis=1).tolist()

          # store reservoirs
          if (ifn==len(atm_exch_files)-mainstream_id):
            reservoirsmask_3d = np.zeros(grid_3d.shape, dtype=bool)
            reservoirsmask_3d[:,:,:] = True
            reservoirsmask_3d_dum = np.zeros(grid_3d.shape, dtype=bool)
            reservoirsmask_3d_dum[:,:,:] = True

            reservoirsmask_3d_dum[np.where(np.logical_and(waterbodyid_grid[:,:,:]>1, waterbodyid_grid[:,:,:]<10000))] = False   
            reservoirsmask_3d[np.where(np.logical_and(reservoirsmask_3d_dum[:,:,:]==False, mask_3d[:,:,:]==False))] = False    
            atm_exch_series["atm_exch_"+specie.get_name().upper()+'_reservoirs']  = np.nansum(np.nansum(np.ma.array(grid_3d_reservoirs, mask=reservoirsmask_3d),axis=2),axis=1).tolist()

        if len(atm_exch_files)>0:
          atm_exch_series["atm_exch_"+specie.get_name().upper()+'_total'] = np.nansum(np.nansum(np.ma.array(tot_spec_atm_exch, mask=mask_3d),axis=2),axis=1).tolist()
      atm_exch_series_list.append(atm_exch_series)
    return atm_exch_series_list

def all_atm_exch_to_table(params):
    atm_exch_series = all_atm_exch_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'atm_exch_'+get_river_name(params)+'.csv')
    dict_to_csv(filename, atm_exch_series)

def all_atm_exch_to_climate_tables(params):
    name_list = ['tropical', 'arid', 'temperate', 'continental', 'polar']
    atm_exch_series_list = atm_exch_to_climate_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))
    i=0
    for atm_exch_series in atm_exch_series_list:
      filename = os.path.join(folder, 'atm_exch_'+name_list[i]+get_river_name(params)+'.csv')
      dict_to_csv(filename, atm_exch_series)
      i+=1

def all_atm_exch_to_pieplots(params):
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'atm_exch_'+get_river_name(params)+'.csv')
    atm_exch_series = csv_to_dict(filename)
    directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))
    emiss_exch_values = np.array([])
    emiss_exch_labels = np.array([])
    uptake_exch_values = np.array([])
    uptake_exch_labels = np.array([])

    for key,value in atm_exch_series.items():
      if (not 'order' in key) and (not key == 'time') and (not 'total' in key):
          val = np.mean(value[int(value.size*0.5):])
          if val > 0:
            emiss_exch_values = np.append(emiss_exch_values, val)
            emiss_exch_labels = np.append(emiss_exch_labels, key.replace("atmospheric_exchange_DIC_", ""))
            emiss_colors = np.append(emiss_colors, global_colors.pop[-1])
          if val < 0:
            uptake_exch_values = np.append(uptake_exch_values, (abs(val)))
            uptake_exch_labels = np.append(uptake_exch_labels, key.replace("atmospheric_exchange_DIC_", ""))
            upt_colors = np.append(upt_colors, global_colors.pop[-1])

    for ival in range(len(emiss_exch_values)):
       emiss_exch_labels[ival] = emiss_exch_labels[ival] + ' ['+'{0:.2f}'.format((emiss_exch_values[ival]/np.sum(emiss_exch_values))*100)+"%]"
    for ival in range(len(uptake_exch_values)):
       uptake_exch_labels[ival] = uptake_exch_labels[ival] + ' ['+'{0:.2f}'.format((uptake_exch_values[ival]/np.sum(uptake_exch_values))*100)+"%]"

    try:
      annotations = ['Results from: \n'+os.path.join(params.outputdir, '..')]
      outname = os.path.join(params.outputdir, '..', 'ANALYSIS', "plots", get_river_name(params)+'_PIECHART_emissions.png')
      make_plots.piechart(emiss_exch_values, emiss_exch_labels, outfilename=outname, title='Carbon emissions from the '+get_river_name(params), annotations=annotations)
    except:
      print('no emissions to plot')

    try:
      annotations = ['Results from: \n'+os.path.join(params.outputdir, '..')]
      outname = os.path.join(params.outputdir, '..', 'ANALYSIS', "plots", get_river_name(params)+'_PIECHART_uptake.png')
      make_plots.piechart(uptake_exch_values, uptake_exch_labels, outfilename=outname, title='Carbon uptake from the '+get_river_name(params), annotations=annotations)    
    except:
      print('no uptake to plot')

def all_atm_exch_to_stackserieplots(params):      
    filename = os.path.join(params.outputdir, '..', 'ANALYSIS', "tables", 'atm_exch_'+get_river_name(params)+'.csv')
    atm_exch_series = csv_to_dict(filename)
    directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))
    time = atm_exch_series['time']
    dellist = list()
    conv = 12*1e-6*(1/params.outputtime)
    for key in atm_exch_series.keys():
        if (('tot' in key) or ('order' in key) or ('time' in key)) and (not 'order7' in key):
            dellist.append(key)
        else:
          atm_exch_series[key] = np.array(atm_exch_series[key])*conv
          atm_exch_series[key][np.isnan(atm_exch_series[key])]=0
          atm_exch_series[key] = atm_exch_series[key].tolist()

    for key in dellist:
            del atm_exch_series[key]

    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'total_budget_'+get_river_name(params)+'.csv')
    budget_series = csv_to_dict(filename)
    y_min, y_max = 0, 1.1*max(budget_series['input'])*conv

    outname = os.path.join(params.outputdir, '..', 'ANALYSIS', "plots", get_river_name(params)+'_TIMESERIES_atm_exch_TgYear.png')	

    make_plots.stacked_linechart(time, atm_exch_series, \
                                           xlab='time [years]', ylab='flux [Tg/year]', \
                                           outfilename=outname, \
                                           title='Atmospheric exchange of CO2 in freshwaters of the '+get_river_name(params)+' basin', \
                                           color='fixed', ylim=(y_min,y_max), xlim=(budget_series['time'][0]-3,budget_series['time'][-1]+3))
def all_atm_exch_per_latitude(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    # read river basin map
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]   

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  

    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]

    waterbodyid.close()
    endo_waterbodyid.close()
    waterbodyoutlet.close()

    lat_array = np.linspace(89.75,-89.75, 360)
    conv = 12*1e-6*(1/params.outputtime)
    ANfolder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "plots"))

    for specie in species:
      atm_exch_files = sorted(directory.get_files_with_str(folder, "atmospheric_exchange_"+specie.get_name().upper()+'*.nc'))
      if len(atm_exch_files):
        sub_spec_atm_exch = np.zeros(np.shape(dummy_nc[dummy_name][-1,:,:]))
        tot_spec_atm_exch = np.zeros(np.shape(dummy_nc[dummy_name][-1,:,:]))
      for ifn in range(len(atm_exch_files)):
        fn = atm_exch_files[ifn]
        nc = Dataset(fn, 'r')
        grid_3d = nc["atmospheric_exchange_"+specie.get_name().upper()][-1,:,:]*conv
        nc.close()

        # store sum of cell
        tot_spec_atm_exch = np.add(grid_3d, tot_spec_atm_exch)

        # store individual orders
        if (ifn < len(atm_exch_files)-mainstream_id) or (ifn > len(atm_exch_files)-mainstream_id):
          atm_lat_array  = np.nansum(np.ma.array(grid_3d, mask=mask_2d),axis=1).tolist()
          filename = os.path.join(ANfolder, 'lat_'+'atm_exch_'+specie.get_val('name')+'_'+str(ifn+1)+'.png')
          make_plots.single_2d(filename, atm_lat_array, lat_array, title="None", ylabel="None", xlabel="None")
          if ifn==6:
            atm_fld_lat_array  = np.nansum(np.ma.array(grid_3d, mask=mask_2d),axis=1).tolist()
        else: # mask the non-outlet lake/reservoirs gridcells in main stream order
          grid_3d_dum = copy.deepcopy(grid_3d)  
          grid_3d_lakes = copy.deepcopy(grid_3d) 
          grid_3d_reservoirs = copy.deepcopy(grid_3d) 
          grid_3d[np.where(waterbodyid_grid[-1,:,:]>1)]=0
          # store main stream without lakes and reservoirs
          atm_rvr_lat_array = np.nansum(np.ma.array(grid_3d, mask=mask_2d),axis=1).tolist()
          filename = os.path.join(ANfolder, 'lat_'+'atm_exch_'+specie.get_val('name')+'_river'+'.png')
          make_plots.single_2d(filename, atm_rvr_lat_array, lat_array, title="None", ylabel="None", xlabel="None")


        # store sum of subgrid orders
        if (ifn < len(atm_exch_files)-mainstream_id):
          sub_spec_atm_exch = np.add(grid_3d, sub_spec_atm_exch)
        elif (ifn == len(atm_exch_files)-mainstream_id):
          atm_sbgrd_lat_array = np.nansum(np.ma.array(sub_spec_atm_exch, mask=mask_2d),axis=1).tolist()
          filename = os.path.join(ANfolder, 'lat_'+'atm_exch_'+specie.get_val('name')+'_subgrid'+'.png')
          make_plots.single_2d(filename, atm_sbgrd_lat_array, lat_array, title="None", ylabel="None", xlabel="None")

        # store lakes
        if (ifn==len(atm_exch_files)-mainstream_id):
          lakesmask_2d = np.zeros(grid_3d.shape, dtype=bool)
          lakesmask_2d[:,:] = True
          lakesmask_2d[np.where(np.logical_and(waterbodyid_grid[-1,:,:]>=10000, mask_2d[:,:]==False))] = False     
          atm_lakes_lat_array = np.nansum(np.ma.array(grid_3d_lakes, mask=lakesmask_2d),axis=1).tolist()
          filename = os.path.join(ANfolder, 'lat_'+'atm_exch_'+specie.get_val('name')+'_lakes'+'.png')
          make_plots.single_2d(filename, atm_lakes_lat_array, lat_array, title="None", ylabel="None", xlabel="None")

        # store reservoirs
        if (ifn==len(atm_exch_files)-mainstream_id):
          reservoirsmask_2d = np.zeros(grid_3d.shape, dtype=bool)
          reservoirsmask_2d[:,:] = True
          reservoirsmask_2d_dum = np.zeros(grid_3d.shape, dtype=bool)
          reservoirsmask_2d_dum[:,:] = True

          reservoirsmask_2d_dum[np.where(np.logical_and(waterbodyid_grid[-1,:,:]>1, waterbodyid_grid[-1,:,:]<10000))] = False   
          reservoirsmask_2d[np.where(np.logical_and(reservoirsmask_2d_dum[:,:]==False, mask_2d[:,:]==False))] = False    
          atm_resrvr_lat_array = np.nansum(np.ma.array(grid_3d_reservoirs, mask=reservoirsmask_2d),axis=1).tolist()    
          filename = os.path.join(ANfolder, 'lat_'+'atm_exch_'+specie.get_val('name')+'_rsrvrs'+'.png')
          make_plots.single_2d(filename, atm_resrvr_lat_array, lat_array, title="None", ylabel="None", xlabel="None")

      if specie.get_val('name')=='DIC':
        atm_tot_lat_array = np.nansum(np.ma.array(tot_spec_atm_exch, mask=mask_2d),axis=1).tolist()    
        filename = os.path.join(ANfolder, 'lat_'+'atm_exch_'+specie.get_val('name')+'_tot'+'.png')
        make_plots.single_2d(filename, atm_tot_lat_array, lat_array, title="None", ylabel="None", xlabel="None")
        filename = os.path.join(ANfolder, 'lat_'+'atm_exch_'+specie.get_val('name')+'_combined'+'.png')
        make_plots.multiple_2d_1frame(filename, [atm_sbgrd_lat_array, atm_rvr_lat_array, atm_lakes_lat_array, atm_resrvr_lat_array, atm_fld_lat_array], [lat_array]*5,\
                           llabel=['small streams', 'large streams', 'lakes', 'reservoirs', 'floodplains/wetlands'], legend=True, xlim="None", ylim="None", title="CO2 exchange from freshwaters per latitude", xlabel="Tg C/yr", ylabel="latitude [degree]", resolution='low')

def all_exports_to_dict(params,add_color=False):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)
    folder = os.path.join(params.outputdir, '..', "BUDGET", "subgrid")

    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")
    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]

    if params.outputtime < 1.: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]-1
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]+1
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
   
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mouthmask_fn = os.path.join(params.water_inputdir, "rivermouth.asc")
    dum_mask = ascraster.create_mask(mouthmask_fn, params.maskid, logical = params.mask_bool_operator, numtype=int)

    if len(dum_mask)>0:
      mouthmask_fn = os.path.join(params.water_inputdir, "rivermouth.asc")
    else: #endoreic basin
      mouthmask_fn = os.path.join(params.water_inputdir, "rivermouths_exporting_to_endoreic_lakes.asc")    

    mouthmask_2d = make_mask.do(mouthmask_fn, params.maskid, dum_asc, mask_type='np_grid',logical=params.mask_bool_operator)
    mouthmask_3d = np.broadcast_to(mouthmask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()

    export_series = dict()
    export_series["time"] = manip.convert_numdate2year(dummy_nc['time'][modeldat_startindex:modeldat_endindex], dummy_nc['time'].units) 

    tot_export = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))

    for specie in species:
      all_export = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
      export_files = directory.get_files_with_str(folder, specie.get_name()+"loadOUT_order6*")
      for fn in export_files:
        nc = Dataset(fn, 'r')
        grid_3d = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
        grid_3d[np.where(nc[specie.get_name()+"loadOUT"][modeldat_startindex:modeldat_endindex,:,:]<1e12)] = nc[specie.get_name()+"loadOUT"][modeldat_startindex:modeldat_endindex,:,:][np.where(nc[specie.get_name()+"loadOUT"][modeldat_startindex:modeldat_endindex,:,:]<1e12)]
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

def all_exports_to_pieplots(params):
    export_series = all_exports_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))
 
    exp_values = np.array([])
    exp_labels = np.array([])
    for key,value_list in export_series.items():
        value = np.array(value_list)
        if (not key == 'time') and (not 'total' in key) and (not 'tss' in key.lower()) and (value.size>0) and (not 'alk' in key.lower()) and (not 'pim' in key.lower()) and (not 'benth' in key.lower()):
          exp_values = np.append(exp_values, np.mean(value[int(value.size*0.5):]))
          exp_labels = np.append(exp_labels, key+' ['+'{0:.2f}'.format((np.mean(value[int(value.size*0.5):])/np.mean(export_series["total_loadOUT"][int(value.size*0.5):]))*100)+"%]")
    annotations = ['Results from: \n'+os.path.join(params.outputdir, '..'), "Average from "+'{0:.2f}'.format(export_series['time'][0])+" to "+'{0:.2f}'.format(export_series['time'][-1])]
    outname = os.path.join(folder, get_river_name(params)+'_PIECHART_export.png')
    make_plots.piechart(exp_values, exp_labels, outfilename=outname, title='Carbon exports from '+get_river_name(params), annotations=annotations, colors='Oranges')

def all_exports_to_stackserieplots(params):   
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))   
    filename = os.path.join(folder, 'exports_'+get_river_name(params)+'.csv')
    export_series = csv_to_dict(filename)
    directory.ensure(os.path.join(folder, '..', "plots"))
    time = export_series['time']
    dellist = list()
    conv = 12*1e-6*(1/params.outputtime)
    for key in export_series.keys():
        if ('ALK' in key) or ('TSS' in key) or ('PIM' in key)or ('time' in key) or ('BENTH' in key) or ('total' in key):
            dellist.append(key)
        else:
            export_series[key] = (np.array(export_series[key])*conv).tolist()
    for key in dellist:
            del export_series[key]

    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'total_budget_'+get_river_name(params)+'.csv')
    budget_series = csv_to_dict(filename)
    y_min, y_max = 0, 1.1*max(budget_series['input'])*conv

    color_dict = dict()
    for key in export_series:
      name = key.replace('_loadOUT')
      color_dict[key] = general_colors[name]
    outname = os.path.join(folder, '..', "plots", get_river_name(params)+'_TIMESERIES_exports.png')	
    make_plots.stacked_linechart(time, export_series, \
                                           xlab='time [years]', ylab='flux [Tg/year]', \
                                           outfilename=outname, \
                                           title='Exports of carbon from '+get_river_name(params),\
                                           color='mixed', color_dict = color_dict, ylim=(y_min, y_max), xlim=(budget_series['time'][0]-3,budget_series['time'][-1]+3))

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
    budget_series['retention'] = ((np.array(budget_series['input'])+np.array(budget_series['emission'])+np.array(budget_series['export']))).tolist()
    return budget_series

def all_budget_to_table(params):
    budget_series = all_budget_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'total_budget_'+get_river_name(params)+'.csv')
    dict_to_csv(filename, budget_series)

def all_budget_to_stackserieplots_Tg_year(params):
    folder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "tables"))
    filename = os.path.join(folder, 'total_budget_'+get_river_name(params)+'.csv')
    budget_series = csv_to_dict(filename)
    line_series = dict()
    line_series['retention'] = (np.array(budget_series['retention'])*12*1e-6*(1/params.outputtime)).tolist()
    time = budget_series['time']
    del budget_series['retention']
    del budget_series['time']
    for key in budget_series:
      budget_series[key]=(np.array(budget_series[key])*12*1e-6*(1/params.outputtime)).tolist()
    y_min, y_max = -1.1*max(budget_series['input']), 1.1*max(budget_series['input'])
    outname = os.path.join(params.outputdir, "..", "ANALYSIS", "plots", get_river_name(params)+'_TIMESERIES_budget.png')	

    make_plots.stacked_linechart_linecombi(time, budget_series, line_series,\
                                           xlab='time [years]', ylab='flux [Tg/yr]', \
                                           outfilename=outname, \
                                           title='Budget of total carbon in '+get_river_name(params) +' basin', \
                                           annotations=[outname], 
                                           color='blue_orange', ylim=(y_min, y_max), xlim=(params.starttime,params.endtime))

def all_fluxes_to_dict(params):
    def list_append(lst, item):
      lst.append(item)
      return lst

    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else: 
      mainstream_id = 1

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  
    

    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]

    mouthmask_fn = os.path.join(params.water_inputdir, "rivermouth.asc")
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mouthmask_2d = make_mask.do(mouthmask_fn, params.maskid, dum_asc, mask_type='np_grid', logical=params.mask_bool_operator)
    mouthmask_3d = np.broadcast_to(mouthmask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    
    flux_series = dict()
    for specie in species:
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
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_smallstreams'] = (np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1)+np.array(flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_subgrid'])).tolist()
              flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_majorstreams'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist() 
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
        flux_series[specie.get_val('name')][specie.get_name()+'_expflux'] = np.nansum(np.nansum(np.ma.array(tot_exp, mask=mouthmask_3d),axis=2),axis=1)  
      flux_series[specie.get_val('name')]['budget'] = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,0,0]))
      for flux in flux_series[specie.get_val('name')]:
        if 'flux' in flux:
          flux_series[specie.get_val('name')]['budget'] += flux_series[specie.get_val('name')][flux]
    dummy_nc.close()
    waterbodyoutlet.close()
    return flux_series

def all_fluxes_per_latitude(params):
    def list_append(lst, item):
      lst.append(item)
      return lst

    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     
    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
    
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1

    lat_array = np.linspace(89.75,-89.75, 360)
    conv = 12*1e-6*(1/params.outputtime)
    ANfolder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "plots"))    

    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]

    mouthmask_fn = os.path.join(params.water_inputdir, "rivermouth.asc")
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    
    for specie in species:
      stoich = reactions.specie_dy(proc,specie.get_name())
      for iproc in range(len(proc)):
        if stoich[iproc]!=0:
          proc_fns = sorted(directory.get_files_with_str(folder, proc[iproc].get_val("name")+"_order*.nc"))
          sub_spec_flux = np.zeros(np.shape(dummy_nc[dummy_name][-1,:,:]))
          tot_spec_flux = np.zeros(np.shape(dummy_nc[dummy_name][-1,:,:]))
          for ifn in range(len(proc_fns)):
            proc_fn = proc_fns[ifn]
            proc_nc = Dataset(proc_fn, 'r')
            grid_3d = proc_nc[proc[iproc].get_val("name")][-1,:,:]*stoich[iproc]*conv
            proc_nc.close()

            # store individual orders
            if (ifn < len(proc_fns)-mainstream_id) or (ifn > len(proc_fns)-mainstream_id):
              order_array = np.nansum(np.ma.array(grid_3d, mask=mask_2d),axis=1).tolist()
              filename = os.path.join(ANfolder, 'lat_'+specie.get_val('name')+'_'+proc[iproc].get_val('name')+'_'+str(ifn)+'.png')
              make_plots.single_2d(filename, order_array, lat_array, title="None", ylabel="None", xlabel="None")
            elif (ifn == len(proc_fns)-mainstream_id): # mask the non-outlet lake/reservoirs gridcells in main stream order
              grid_3d_dum = copy.deepcopy(grid_3d) 
              grid_3d_dum[np.where(waterbodyid_grid[-1,:,:]>1)]=0
			  # store main stream without lakes and reservoirs
              river_array = np.nansum(np.ma.array(grid_3d_dum, mask=mask_2d),axis=1).tolist()
              filename = os.path.join(ANfolder, 'lat_'+specie.get_val('name')+'_'+proc[iproc].get_val('name')+'_river'+'.png')
              make_plots.single_2d(filename, river_array, lat_array, title="None", ylabel="None", xlabel="None")

            # store sum of subgrid orders
            if (ifn < len(proc_fns)-mainstream_id):
              sub_spec_flux = np.add(grid_3d, sub_spec_flux)
            elif (ifn == len(proc_fns)-mainstream_id):
               sbgrd_array = np.nansum(np.ma.array(sub_spec_flux, mask=mask_2d),axis=1).tolist()
               filename = os.path.join(ANfolder, 'lat_'+specie.get_val('name')+'_'+proc[iproc].get_val('name')+'_sbgrd'+'.png')
               make_plots.single_2d(filename, sbgrd_array, lat_array, title="None", ylabel="None", xlabel="None")
    
            # store sum of cell
            tot_spec_flux = np.add(grid_3d, tot_spec_flux)

            # store lakes
            if (ifn==len(proc_fns)-mainstream_id):
              lakesmask_2d = np.zeros(grid_3d.shape, dtype=bool)
              lakesmask_2d[:,:] = True
              lakesmask_2d[np.where(np.logical_and(waterbodyid_grid[-1,:,:]>=10000, mask_2d[:,:]==False))] = False     
              lakes_array = np.nansum(np.ma.array(grid_3d, mask=lakesmask_2d),axis=1).tolist()
              filename = os.path.join(ANfolder, 'lat_'+specie.get_val('name')+'_'+proc[iproc].get_val('name')+'_lakes'+'.png')
              make_plots.single_2d(filename, lakes_array, lat_array, title="None", ylabel="None", xlabel="None")

            # store reservoirs
            if (ifn==len(proc_fns)-mainstream_id):
              reservoirsmask_2d = np.zeros(grid_3d.shape, dtype=bool)
              reservoirsmask_2d[:,:] = True
              reservoirsmask_2d_dum = np.zeros(grid_3d.shape, dtype=bool)
              reservoirsmask_2d_dum[:,:] = True

              reservoirsmask_2d_dum[np.where(np.logical_and(waterbodyid_grid[-1,:,:]>1, waterbodyid_grid[-1,:,:]<10000))] = False   
              reservoirsmask_2d[np.where(np.logical_and(reservoirsmask_2d_dum[:,:]==False, mask_2d[:,:]==False))] = False   
              reservoir_array = np.nansum(np.ma.array(grid_3d, mask=mask_2d),axis=1).tolist()
              filename = os.path.join(ANfolder, 'lat_'+specie.get_val('name')+'_'+proc[iproc].get_val('name')+'_rsrvrs'+'.png')
              make_plots.single_2d(filename, reservoir_array, lat_array, title="None", ylabel="None", xlabel="None")

            # store floodplains
            if 'order7' in proc_fns[ifn]:
              fldpln_array = np.nansum(np.ma.array(grid_3d, mask=mask_2d),axis=1).tolist()     
              filename = os.path.join(ANfolder, 'lat_'+specie.get_val('name')+'_'+proc[iproc].get_val('name')+'_fldpln'+'.png')
              make_plots.single_2d(filename, fldpln_array, lat_array, title="None", ylabel="None", xlabel="None")
  
          tot_array = np.nansum(np.ma.array(tot_spec_flux, mask=mask_2d),axis=1).tolist() 
          filename = os.path.join(ANfolder, 'lat_'+specie.get_val('name')+'_'+proc[iproc].get_val('name')+'_tot'+'.png')
          make_plots.single_2d(filename, tot_array, lat_array, title="None", ylabel="None", xlabel="None")

          filename = os.path.join(ANfolder, 'lat_'+specie.get_val('name')+'_'+proc[iproc].get_val('name')+'_combined'+'.png')
          make_plots.multiple_2d_1frame(filename, [sbgrd_array, river_array, lakes_array, reservoir_array, fldpln_array], [lat_array]*5,\
                           llabel=['small streams', 'large streams', 'lakes', 'reservoirs', 'floodplains/wetlands'], legend=True, xlim="None", ylim="None", title=proc[iproc].get_val('name')+" from freshwaters per latitude", xlabel="Tg C/yr", ylabel="latitude [degree]", resolution='low')


def all_fluxes_to_climate_dict(params):
    def list_append(lst, item):
      lst.append(item)
      return lst

    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index_species.make_index_species(params,species,proc)

    if params.lfloodplains:
      mainstream_id = 2
    else:
      mainstream_id = 1

    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

    dummy_nc = Dataset(proclist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     
    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d_dum = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d_dum = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  
    
    print('modeldat_startindex=', modeldat_startindex)
    print('modeldat_endindex=', modeldat_endindex)
    print('all_dat_startindex=', all_dat_startindex)
    print('modeldat_endindex=', modeldat_endindex)

    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]

    mouthmask_fn = os.path.join(params.water_inputdir, "rivermouth.asc")
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mouthmask_2d = make_mask.do(mouthmask_fn, params.maskid, dum_asc, mask_type='np_grid', logical=params.mask_bool_operator)
    mouthmask_3d = np.broadcast_to(mouthmask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()

    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    tropics_mask_2d = make_mask.do(climatemask_fn, 1, dum_asc, mask_type='np_grid', logical='EQ')
    tropics_mask_3d = np.broadcast_to(tropics_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    arid_mask_2d = make_mask.do(climatemask_fn, 2, dum_asc, mask_type='np_grid', logical='EQ')
    arid_mask_3d = np.broadcast_to(arid_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    temperate_mask_2d = make_mask.do(climatemask_fn, 3, dum_asc, mask_type='np_grid', logical='EQ')
    temperate_mask_3d = np.broadcast_to(temperate_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    continental_mask_2d = make_mask.do(climatemask_fn, 4, dum_asc, mask_type='np_grid', logical='EQ')
    continental_mask_3d = np.broadcast_to(continental_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    polar_mask_2d = make_mask.do(climatemask_fn, 5, dum_asc, mask_type='np_grid', logical='EQ')
    polar_mask_3d = np.broadcast_to(polar_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
 
    name_list = ['tropical', 'arid', 'temperate', 'continental', 'polar']
    climate_mask_list = [tropics_mask_3d, arid_mask_3d, temperate_mask_3d, continental_mask_3d, polar_mask_3d]
    flux_series_list = []
    for climate_name in name_list:
      mask_3d = np.zeros(mask_3d_dum.shape, dtype=bool)
      mask_3d[:,:,:] = True
      mask_3d[np.logical_and(climate_mask_list[name_list.index(climate_name)][:,:,:]==False, mask_3d_dum[:,:,:]==False)] = False
      flux_series = dict()
      for specie in species:
       if 'dic' in specie.get_val('name').lower():
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
                flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_smallstreams'] = (np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1)+np.array(flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_subgrid'])).tolist()
                flux_series[specie.get_val('name')][proc[iproc].get_val("name")+'_majorstreams'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist() 
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
          flux_series[specie.get_val('name')][specie.get_name()+'_expflux'] = np.nansum(np.nansum(np.ma.array(tot_exp, mask=mouthmask_3d),axis=2),axis=1)  
        flux_series[specie.get_val('name')]['budget'] = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,0,0]))
        #for flux in flux_series[specie.get_val('name')]:
        #  if 'flux' in flux:
        #    flux_series[specie.get_val('name')]['budget'] += flux_series[specie.get_val('name')][flux]
      flux_series_list.append(flux_series)
    dummy_nc.close()
    waterbodyoutlet.close()
    return flux_series_list

def all_fluxes_to_table(params):
    flux_series = all_fluxes_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))
    for speciekey, speciedict in flux_series.items():
      filename = os.path.join(folder, speciekey+'_biogeochemistry_'+get_river_name(params)+'.csv')
      dict_to_csv(filename, speciedict)

def all_fluxes_to_climate_tables(params):
    name_list = ['tropical', 'arid', 'temperate', 'continental', 'polar']
    flux_series_list = all_fluxes_to_climate_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))
    i=0
    for flux_series in flux_series_list:
      for speciekey, speciedict in flux_series.items():
        filename = os.path.join(folder, speciekey+'_'+name_list[i]+'_biogeochemistry_'+get_river_name(params)+'.csv')
        dict_to_csv(filename, speciedict)
      i+=1

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

    filelist = directory.get_files_with_str(folder, arguments[0]+"*_order6*")

    dummy_nc = Dataset(filelist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(filelist[0]))[0][:-7]

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     
    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  

    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]    

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
                sec_series[arg][arg+'_smallstreams'] = (np.add(np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_small_streams_3d),axis=2),axis=1),np.array(sec_series[arg][arg+'_subgrid']))).tolist()
                sec_series[arg][arg+'_majorstreams'] = np.nansum(np.nansum(np.ma.array(grid_3d, mask=mask_major_streams_3d),axis=2),axis=1).tolist() 
                grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]
     
            # store sum of cell
            grid_3d[np.where(np.isnan(grid_3d))] = 0
            tot_arg = np.add(grid_3d,tot_arg)

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
    waterbodyoutlet.close()
    return sec_series

def all_sec_to_climate_dict(params):
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

    filelist = directory.get_files_with_str(folder, arguments[0]+"*_order6*")

    dummy_nc = Dataset(filelist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(filelist[0]))[0][:-7]

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     
    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d_dum = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d_dum = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  

    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]   
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    tropics_mask_2d = make_mask.do(climatemask_fn, 1, dum_asc, mask_type='np_grid', logical='EQ')
    tropics_mask_3d = np.broadcast_to(tropics_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    arid_mask_2d = make_mask.do(climatemask_fn, 2, dum_asc, mask_type='np_grid', logical='EQ')
    arid_mask_3d = np.broadcast_to(arid_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    temperate_mask_2d = make_mask.do(climatemask_fn, 3, dum_asc, mask_type='np_grid', logical='EQ')
    temperate_mask_3d = np.broadcast_to(temperate_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    continental_mask_2d = make_mask.do(climatemask_fn, 4, dum_asc, mask_type='np_grid', logical='EQ')
    continental_mask_3d = np.broadcast_to(continental_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
    polar_mask_2d = make_mask.do(climatemask_fn, 5, dum_asc, mask_type='np_grid', logical='EQ')
    polar_mask_3d = np.broadcast_to(polar_mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape).copy()
 
    name_list = ['tropical', 'arid', 'temperate', 'continental', 'polar']
    climate_mask_list = [tropics_mask_3d, arid_mask_3d, temperate_mask_3d, continental_mask_3d, polar_mask_3d]
    sec_series_list = []
    for climate_name in name_list:
      mask_3d = np.zeros(mask_3d_dum.shape, dtype=bool)
      mask_3d[:,:,:] = True
      mask_3d[np.logical_and(climate_mask_list[name_list.index(climate_name)][:,:,:]==False, mask_3d_dum[:,:,:]==False)] = False
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
      sec_series_list.append(sec_series)
    return sec_series_list

def all_sec_per_latitude(params):
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

    filelist = directory.get_files_with_str(folder, arguments[0]+"*_order6*")

    dummy_nc = Dataset(filelist[0], 'r')
    dummy_name = os.path.splitext(os.path.basename(filelist[0]))[0][:-7]

    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     

    if params.outputtime < 1: 
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101_mon.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101_mon.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101_mon.nc"), 'r')  
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)     
    else:
      waterbodyoutlet  = Dataset(os.path.join(params.water_inputdir, "waterbodyoutlet_101.nc"), 'r')
      waterbodyid = Dataset(os.path.join(params.water_inputdir, "waterbodyid_101.nc"), 'r')
      endo_waterbodyid = Dataset(os.path.join(params.water_inputdir, "endo_waterbodyid_101.nc"), 'r')
      modelrun_start = max(dummy_nc['time'][0],0)
      modelrun_end = dummy_nc['time'][-1]
      all_dat_startindex = np.where(waterbodyoutlet['time'][:] >= modelrun_start)[0][0]
      all_dat_startindex=max(all_dat_startindex, 0)
      all_dat_endindex = np.where(waterbodyoutlet['time'][:] <= modelrun_end)[0][-1]
      modeldat_startindex = np.where(dummy_nc['time'][:] >= waterbodyoutlet['time'][all_dat_startindex])[0][0]
      modeldat_endindex = np.where(dummy_nc['time'][:] <= waterbodyoutlet['time'][all_dat_endindex])[0][-1]+1
      if modeldat_endindex==len(dummy_nc['time'][:]):
        all_dat_endindex +=1
      mask_3d = np.broadcast_to(mask_2d, dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:].shape)  

    waterbodyoutlet_grid = waterbodyoutlet['waterbodyoutlet'][all_dat_startindex:all_dat_endindex,:,:]
    waterbodyid_grid = waterbodyid['waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]
    endo_waterbodyid_grid = endo_waterbodyid['endo_waterbodyid'][all_dat_startindex:all_dat_endindex,:,:]    

    sec_series = dict()
    for arg in arguments:
      sec_series[arg] = dict()
      sec_series[arg]["time"] = manip.convert_numdate2year(dummy_nc['time'][:], dummy_nc['time'].units)[modeldat_startindex:modeldat_endindex]
      arg_fns = sorted(directory.get_files_with_str(folder, arg+"_order*.nc"))
      sub_arg = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
      tot_arg = np.zeros(np.shape(dummy_nc[dummy_name][modeldat_startindex:modeldat_endindex,:,:]))
      for ifn in range(len(arg_fns)):
            arg_fn = arg_fns[ifn]
            arg_nc = Dataset(arg_fn, 'r')
            grid_3d = arg_nc[arg][modeldat_startindex:modeldat_endindex,:,:]
            # store individual orders
            if (ifn < len(arg_fns)-mainstream_id) or (ifn > len(arg_fns)-mainstream_id):
              sec_series[arg][arg+'_order'+str(arg_fn[-4])] = np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2).tolist()
            elif (ifn == len(arg_fns)-mainstream_id): # mask the non-outlet lake/reservoirs gridcells in main stream order
              grid_3d_dum = copy.deepcopy(grid_3d) 
              grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]=0
			  # store main stream without lakes and reservoirs
              sec_series[arg][arg+'_river'] = np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2).tolist()
			  
              grid_3d[np.where(waterbodyid_grid[:,:,:]>1)]+=grid_3d_dum[np.where(waterbodyid_grid[:,:,:]>1)]
              sec_series[arg][arg+'_order'+str(arg_fn[-4])] = np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2).tolist()

            
            # store sum of subgrid orders
            if (ifn < len(arg_fns)-mainstream_id):
              sub_arg = np.add(sub_arg, np.ma.array(grid_3d, mask=mask_3d)*params.number_of_rivers[ifn])
            elif (ifn == len(arg_fns)-mainstream_id):
              sec_series[arg][arg+'_subgrid'] = np.nansum(np.ma.array(sub_arg, mask=mask_3d),axis=2).tolist()
     
            # store sum of cell
            tot_arg = np.add(grid_3d, tot_arg)

            # store lakes
            if (ifn==len(arg_fns)-mainstream_id):
              lakesmask_3d = np.zeros(grid_3d.shape, dtype=bool)
              lakesmask_3d[:,:,:] = True
              lakesmask_3d[np.where(np.logical_and(waterbodyid_grid[:,:,:]>=10000, mask_3d[:,:,:]==False))] = False     
              sec_series[arg][arg+'_lakes']  = np.nansum(np.ma.array(grid_3d, mask=lakesmask_3d),axis=2).tolist()

            # store reservoirs
            if (ifn==len(arg_fns)-mainstream_id):
              reservoirsmask_3d = np.zeros(grid_3d.shape, dtype=bool)
              reservoirsmask_3d[:,:,:] = True
              reservoirsmask_3d_dum = np.zeros(grid_3d.shape, dtype=bool)
              reservoirsmask_3d_dum[:,:,:] = True

              reservoirsmask_3d_dum[np.where(np.logical_and(waterbodyid_grid[:,:,:]>1, waterbodyid_grid[:,:,:]<10000))] = False   
              reservoirsmask_3d[np.where(np.logical_and(reservoirsmask_3d_dum[:,:,:]==False, mask_3d[:,:,:]==False))] = False   
              sec_series[arg][arg+'_reservoirs']  = np.nansum(np.ma.array(grid_3d, mask=reservoirsmask_3d),axis=2).tolist()

            # store floodplains
            if 'order7' in arg_fns[ifn]:
              sec_series[arg][arg+'_floodplains/wetlands']  = np.nansum(np.ma.array(grid_3d, mask=mask_3d),axis=2).tolist() 
              sec_series[arg].pop(arg+'_order7', None)             

  
      sec_series[arg][arg+'_tot'] = np.nansum(np.ma.array(tot_arg, mask=mask_3d),axis=2).tolist()
    dummy_nc.close()
    waterbodyoutlet.close()

    lat_array = np.linspace(89.75,-89.75, 360)
    ANfolder = directory.ensure(os.path.join(params.outputdir, "..", "ANALYSIS", "plots"))    


    for arg in sec_series.keys():
      if arg=='vol':
        unit='km3'
        xlim = (0, 400)
      elif arg=='area':
        unit='km2'
        xlim = (0, 80000)
      for itime in range(len(sec_series[arg]['time'])):
        sbgrd_array = sec_series[arg][arg+'_subgrid'][itime]
        river_array = sec_series[arg][arg+'_river'][itime]
        lakes_array = sec_series[arg][arg+'_lakes'][itime]
        reservoir_array = sec_series[arg][arg+'_reservoirs'][itime]
        fldpln_array = sec_series[arg][arg+'_floodplains/wetlands'][itime]

        filename = os.path.join(ANfolder, 'lat_'+arg+'_combined_'+str(sec_series[arg]['time'][itime])+'.png')
        make_plots.multiple_2d_1frame(filename, [sbgrd_array, river_array, lakes_array, reservoir_array, fldpln_array], [lat_array]*5,\
                           llabel=['small streams', 'large streams', 'lakes', 'reservoirs', 'floodplains/wetlands'], legend=True, xlim=xlim, ylim="None", title=arg+" of freshwaters per latitude", xlabel=unit, ylabel="latitude [degree]", resolution='low')

      outfilename = os.path.join(ANfolder, 'lat_'+arg+'_combined'+'.gif') 
      images=list()
      filenames = sorted(directory.get_files_with_str(ANfolder, 'lat_'+arg+"*.png"))
      for filename in filenames:
        images.append(imageio.imread(filename))
      imageio.mimsave(outfilename, images, duration=0.3)

def all_stream_env_conditions_to_climate_tables(params):
    name_list = ['tropical', 'arid', 'temperate', 'continental', 'polar']
    sec_series_list = all_sec_to_climate_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))
    i=0
    for sec_series in sec_series_list:
      for key, dictionary in sec_series.items():
        filename = os.path.join(folder, key+'_'+name_list[i]+'_sec_'+get_river_name(params)+'.csv')
        dict_to_csv(filename, dictionary)
      i+=1

def all_stream_env_conditions_to_table(params):
    sec_series = all_sec_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))
    for key, dictionary in sec_series.items():
      filename = os.path.join(folder, key+'_sec_'+get_river_name(params)+'.csv')
      dict_to_csv(filename, dictionary)

def all_fluxes_to_multiplots(params):
    try:
      species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
      flux_series = dict()
      for specie in species:
        filename = os.path.join(params.outputdir, '..', 'ANALYSIS', "tables", specie.get_val('name')+'_biogeochemistry_'+get_river_name(params)+'.csv')
        flux_series[specie.get_val('name')] = csv_to_dict(filename)
    except:
      flux_series = all_fluxes_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))

    multiplotdict = dict()
    for speciename, speciedict in flux_series.items():
      for key, array in speciedict.items():
        if ('subgrid' in key) or ('river' in key) or ('lakes' in key) or ('reservoirs' in key) or ('totflux' in key) or ('order7' in key):
          key = key.replace('order7', 'floodplains')
          procname = key.replace('_subgrid','').replace('_river', '').replace('_lakes', '').replace('_reservoirs', '').replace('_totflux', '').replace('_floodplains', '')
          print(key)
          if (procname not in multiplotdict):
              multiplotdict[procname] = dict()
          if (key.replace(procname+'_', '') not in multiplotdict[procname]):
              multiplotdict[procname][key.replace(procname+'_', '')] = [abs(number) for number in array]
    for procname, componentsdict in multiplotdict.items():
        lx = list()
        ly = list()
        ltitles = list()
        outfilename = os.path.join(folder, procname+'_'+get_river_name(params)+'_timeseries.png')
        for component, array in componentsdict.items():
          lx.append(speciedict['time'])
          ly.append(array)
          ltitles.append(component)
        make_plots.multiframe_1fig(outfilename, lx, ly, suptitle=procname+' in the '+get_river_name(params)+' basin', ltitles=ltitles, lunits='None', lxlabels='None', linestyle = '-', marker=' ', markersize = 1)

def all_speciefluxes_to_multiplots(params):
    try:
      species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
      flux_series = dict()
      for specie in species:
        filename = os.path.join(params.outputdir, '..', 'ANALYSIS', "tables", specie.get_val('name')+'_biogeochemistry_'+get_river_name(params)+'.csv')
        flux_series[specie.get_val('name')] = csv_to_dict(filename)
    except:
      flux_series = all_fluxes_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS',"plots"))

    multiplotdict = dict()
    for speciename, speciedict in flux_series.items():
      lx = list()
      ly = list()
      ltitles = list()
      outfilename = os.path.join(folder, speciename+'_'+get_river_name(params)+'_timeseries.png')
      for key, array in speciedict.items():
        if ('totflux' in key):
          procname = key.replace('_totflux', '')
          lx.append(speciedict['time'])
          ly.append(array)
          ltitles.append(procname)
      if len(lx)>0:
        make_plots.multiframe_1fig(outfilename, lx, ly, suptitle=speciename+' in the '+get_river_name(params)+' basin', ltitles=ltitles, lunits='None', lxlabels='None', linestyle = '-', marker=' ', markersize = 1)

def all_fluxes_to_pieplots(params):
    try:
      species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
      flux_series = dict()
      for specie in species:
        filename = os.path.join(params.outputdir, '..', 'ANALYSIS', "tables", specie.get_val('name')+'_biogeochemistry_'+get_river_name(params)+'.csv')
        flux_series[specie.get_val('name')] = csv_to_dict(filename)
    except:
      flux_series = all_fluxes_to_dict(params)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))

    title_list = ['origins','fates']
    for speciekey, speciedict in flux_series.items():
      pos_list = []
      neg_list = []
      pos_labels = []
      neg_labels = [] 
      for key, values in speciedict.items():
        if ('flux' in key):
          if np.mean(values)>0:
            pos_list.append(np.mean(values))
            pos_labels.append(key.replace('flux', ''))
          elif np.mean(values)<0:
            neg_list.append(abs(np.mean(values)))
            neg_labels.append(key.replace('flux', ''))
			
      data_list = [pos_list, neg_list]
      labels = [pos_labels, neg_labels]
      if len(data_list[1])>0 or len(data_list[0])>0:
        outname = os.path.join(folder, get_river_name(params)+'_PIECHARTS_'+speciekey+'_origins_fates.png')		
        make_plots.multi_piechart(data_list, labels, title_list, outfilename=outname, title=speciekey + ' in ' + get_river_name(params) + ' freshwaters ', annotations=["no comments"])

def all_fluxes_per_wbtype_to_pieplots(params):
    #try:
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    flux_series = dict()
    for specie in species:
      filename = os.path.join(params.outputdir, '..', 'ANALYSIS', "tables", specie.get_val('name')+'_biogeochemistry_'+get_river_name(params)+'_Tg.csv')
      print(filename)
      flux_series[specie.get_val('name')] = csv_to_dict(filename)
    #except:
    #  flux_series = all_fluxes_to_dict(params)

    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))

    for speciekey, speciedict in flux_series.items():
        pie_dict = {}    
        for key, values in speciedict.items():

          if ('lakes' in key) or ('reservoirs' in key) or ('subgrid' in key) or ('river' in key) or ('order7' in key):
            procname = '_'.join(key.split('_')[:-1])
            if not procname in pie_dict.keys():
              pie_dict[procname] = {}
              pie_dict[procname]['labels'] = []
              pie_dict[procname]['data'] = []
              pie_dict[procname]['labels'].append(key.split('_')[-1].replace('order7', 'floodplain'))
              pie_dict[procname]['data'].append(abs(sum(values)/len(values)))
            else:
              pie_dict[procname]['labels'].append(key.split('_')[-1].replace('order7', 'floodplain'))
              try:
                pie_dict[procname]['data'].append(abs(sum(values)/len(values)))
              except:
                print("sum(values)=", sum(values))
                print("len(values)=", len(values))
        for procname in pie_dict.keys():
          outname = os.path.join(folder, get_river_name(params)+'_PIECHARTS_perwaterbody_'+procname+'.png')
          try:		
            make_plots.piechart(pie_dict[procname]['data'], pie_dict[procname]['labels'], outfilename=outname, title=procname+'per waterbody type', annotations=["no comments"], colors='Oranges')   
          except(ValueError):
            pass

def all_fluxes_per_basin_to_pieplots(params):
    pass

def all_fluxes_per_climate_to_pieplots(params):
    pass

def all_fluxes_to_stackserieplots(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    flux_series = dict()
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "tables"))    
    for specie in species:
        filename = os.path.join(folder, specie.get_val('name')+'_biogeochemistry_'+get_river_name(params)+'.csv')
        flux_series[specie.get_val('name')] = csv_to_dict(filename)
    #flux_series = all_fluxes_to_dict(params)
    directory.ensure(os.path.join(folder, '..', "plots"))

    filename = os.path.join(folder, 'total_budget_'+get_river_name(params)+'.csv')
    budget_series = csv_to_dict(filename)
    y_min, y_max = -1.1*max(budget_series['srcflux']), 1.1*max(budget_series['srcflux'])
    
    for speciekey, speciedict in flux_series.items():
        if len(speciedict.keys())>0:
          time = speciedict['time']
          dellist = list()
          for key in speciedict.keys():
            if (not 'totflux' in key) and (not 'expflux' in key) and (not 'src' in key):
              dellist.append(key)
          for key in dellist:
            del speciedict[key]

          outname = os.path.join(folder, '..', "plots", get_river_name(params)+'_Mmol_TIMESERIES_'+speciekey+'_origins_fates.png')	
          make_plots.stacked_linechart(time, speciedict, \
                                           xlab='time [years]', ylab='flux [Mmol/month]', \
                                           outfilename=outname, \
                                           title='Biogeochemical and lateral fluxes of '+speciekey+' in the '+ get_river_name(params), \
                                           annotations=[outname], 
                                           color='blue_orange', ylim=(y_min, y_max), xlim=(params.starttime,params.endtime))

def all_fluxes_to_stackserieplots_TgPerYear(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    folder = directory.ensure(os.path.join(params.outputdir, '..', 'ANALYSIS', "plots"))
    flux_series = dict()
    line_series = dict()
    for specie in species:
        filename = os.path.join(folder, '..', "tables", specie.get_val('name')+'_biogeochemistry_'+get_river_name(params)+'.csv')
        flux_series[specie.get_val('name')] = csv_to_dict(filename)
    #flux_series = all_fluxes_to_dict(params)

    filename = os.path.join(folder, '..', "tables", 'total_budget_'+get_river_name(params)+'.csv')
    budget_series = csv_to_dict(filename)
    y_min, y_max = -1.1*max(budget_series['input'])*12*1e-6*(1/params.outputtime), 1.1*max(budget_series['input'])*12*1e-6*(1/params.outputtime)

    for speciekey, speciedict in flux_series.items():
        if len(speciedict.keys())>0:
          time = speciedict['time']

          dellist = list()
          for key in speciedict.keys():
            if (not 'totflux' in key) and (not 'expflux' in key) and (not 'src' in key) and (not 'budget' in key):
              dellist.append(key)
          for key in dellist:
            del speciedict[key]
          speciedict2 = dict()
          for key, array in speciedict.items():
            speciedict2[key] = (np.array(array)*(1e-6)*(1/params.outputtime)*12).tolist()
          line_series['budget'] = (np.array(speciedict2['budget'])).tolist()
          del speciedict2['budget']

          outname = os.path.join(folder, get_river_name(params)+'_Tg_TIMESERIES_'+speciekey+'_origins_fates.png')	
          make_plots.stacked_linechart_linecombi(time, speciedict2, line_series,\
                                           xlab='time [years]', ylab='flux [Tg/yr]', \
                                           outfilename=outname, \
                                           title='Budget of '+ speciekey + ' in the '+get_river_name(params) +' basin', \
                                           annotations=[outname], 
                                           color='blue_orange', ylim=(y_min, y_max), xlim=(params.starttime, params.endtime))

def piechart_export(params, loc=[None]):
    print('piechart_export start')
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index.species(params,species,proc)
    
    basin = ascraster.Asciigrid(ascii_file=params.basin,numtype=int)
    icell = manip.calc_index_from_coordin(loc[0], loc[1], basin.nrows, basin.ncols, basin.xllcorner, basin.yllcorner, basin.cellsize)
    ilat, ilon = manip.calc_row_col_from_index(icell, basin.ncols)
    discharge = Dataset(os.path.join(params.water_inputdir,params.discharge), 'r')[params.discharge_varname][0,ilat,ilon]
    
    labels = list()
    carbon_array = np.array([])
    mass_transport = np.array([])
    cols=list()
    for specie in species: #for every specie
        if 'alk' not in specie.get_name().lower() and 'benth' not in specie.get_name():
            searchdir = params.states_outputdir
            flist = directory.get_files_with_str(searchdir, "conc_"+specie.get_name()+'.nc')
            labels.append(specie.get_name())
            for f in flist:

                nc_dat = Dataset(f, 'r')
                var = "conc_"+specie.get_name()
                carbon_array = np.append(carbon_array, nc_dat[var][0,ilat,ilon])

                #teragrams per year
                mass_transport = np.append(mass_transport, carbon_array[-1]*discharge*1e-3)
                cols.append(specie.get_val('color'))

    colors = [x for _,x in sorted(zip(mass_transport,cols), reverse=True)]

    for i in range(len(carbon_array)): #for every specie
        labels[i]=labels[i]+' '+'{0:.4f}'.format(mass_transport[i])+' Tg/yr  ('+str(int(carbon_array[i]*100/carbon_array.sum()))+'%)'
        
    annotations = ['Results from: \n'+os.path.join(params.outputdir, '..')]
    outname = os.path.join(params.outputdir, '..', 'ANALYSIS',  str(params.maskid)+'_export.png')
    make_plots.piechart(mass_transport, labels, outfilename=outname, title='Carbon export of the '+str(params.maskid), annotations=annotations, colors=colors, xlim=(params.starttime,params.endtime))   
    print('piechart_export finish')

def map_net_budget(params):
    species,sources,proc,params_local = read_parameter.readfile(params.species_ini)
    make_index.species(params,species,proc)
    basinid = pickle.load(open(os.path.join(params.outputdir,'mask_asc.pkl'),'rb'))
    folder = os.path.join(params.outputdir, "..", "BUDGET", "subgrid")
    dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
    mask_2d_dum = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')
    climatemask_fn = os.path.join(params.water_inputdir, 'climate.asc')
    climate_mask_2d_dum = make_mask.do(climatemask_fn, 0, dum_asc, logical='GT', mask_type='np_grid')
    mask_2d = np.zeros(mask_2d_dum.shape, dtype=bool)
    mask_2d[:,:] = True
    mask_2d[np.where(np.logical_and(mask_2d_dum[:,:]==False, climate_mask_2d_dum[:,:]==False))] = False     
    try:
        all_budget
    except NameError:
        proclist = directory.get_files_with_str(folder, species[0].get_name().upper()+"*_order6*")

        dummy_nc = Dataset(proclist[0], 'r')
        dummy_name = os.path.splitext(os.path.basename(proclist[0]))[0][:-7]    

        specie_budget = np.zeros(np.shape(dummy_nc[dummy_name]))
        all_budget = np.zeros(np.shape(dummy_nc[dummy_name]))
    else:
        pass

    for specie in species:
        stoich = reactions.specie_dy(proc,specie.get_name())
        specie_budget = np.zeros(np.shape(dummy_nc[dummy_name]))
        for iproc in range(len(proc)):
            proc_fn = os.path.join(folder, proc[iproc].get_val("name")+"_order6.nc")
            proc_nc = Dataset(proc_fn, 'r')
            grid_3d = proc_nc[proc[iproc].get_val("name")][:,:,:]*stoich[iproc]
            specie_budget = np.add(specie_budget, grid_3d)
        # also account for advection terms
        adv_terms = ['upstrloadIN','loadOUT','srcloadIN']
        adv_pos = [1,-1, 1]
        for i in range(len(adv_terms)):
            proc_fn = os.path.join(folder, specie.get_name()+'_'+adv_terms[i]+"_order6.nc")
            proc_nc = Dataset(proc_fn, 'r')
            grid_3d = proc_nc[specie.get_name()+'_'+adv_terms[i]][:,:,:]*adv_pos[i]
            specie_budget = np.add(specie_budget, grid_3d)
        mask_3d = np.broadcast_to(mask_2d, specie_budget.shape)
        specie_budget = np.ma.array(specie_budget, mask=mask_3d)
        if ((specie.get_name().lower()!='alk') and (not 'tss' in specie.get_name().lower())):
          all_budget += np.add(all_budget, specie_budget)
        
    
        specie_budget_nc = Dataset(os.path.join(folder, specie.get_name()+"_budget.nc"), 'w')
        output_conversion.init_ncdata(folder, specie_budget_nc, specie.get_name()+"_budget", basinid, unit='Mmol/yr', long_name=specie.get_name()+"_budget")
        for itime in range(len(dummy_nc['time'])):
            manip.add_grid_time(specie_budget_nc, specie.get_name()+"_budget", specie_budget[itime,:,:], dummy_nc['time'][itime])

    all_budget = np.ma.array(all_budget, mask=mask_3d)
    all_budget_nc = Dataset(os.path.join(folder, "all_budget.nc"), 'w')
    output_conversion.init_ncdata(folder, all_budget_nc, "all_budget", basinid, unit='Mmol/yr', long_name="all_budget")
    for itime in range(len(dummy_nc['time'])):
        manip.add_grid_time(all_budget_nc, "all_budget", all_budget[itime,:,:], dummy_nc['time'][itime])

if __name__ == "__main__":
    # Set the general path for the own python modules

    import general_path   
    print(sys.path)
    # Parse command-line arguments and set parameters for script
    #try:
    # Parse command-line arguments and set parameters for script
    # Startup logging and runtime start
    params,log,s = general_startup.general_startup(sys.argv)
    #except SystemExit:
    #  raise MyError("Error has occured in the reading of the commandline options.")
    
    do(params)	
