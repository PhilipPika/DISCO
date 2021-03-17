from datetime import datetime
import os
from netCDF4 import Dataset
import numpy as np
import sys

import extract_varname
import make_mask
print('HELLO HELO')
import general_path
import ascraster
import output_conversion
import manip

def make_time_array(length, unit='hours since 1900-01-01 00:00:00', startyear=1900, startmonth=1, startday=1, timestep='year'):
    '''
    function to make a time array with:
       length 'length' as int
       unit 'unit' as defined in netCDF file
       starting year 'startyear' as int
       starting month 'startmonth' as int
       starting day 'startday' as int
       timestep 'timestep as string. OPTIONS: 'year', 'month'
    '''
    time_array=[]

    if timestep=='year':
      for itime in range(length):
        timestr = str(startyear+itime)+'-'+str(startmonth)+'-01T00:00:00Z'
        dt_time = datetime.strptime(timestr, '%Y-%m-%dT%H:%M:%SZ')
        if unit=='hours since 1900-01-01 00:00:00':
          ini_timestr = str(1900)+'-'+str(1)+'-'+str(1)+'T00:00:00Z'
          ini_dt_time = datetime.strptime(ini_timestr, '%Y-%m-%dT%H:%M:%SZ')
          td = dt_time - ini_dt_time
          time_array.append(td.total_seconds()//3600)
    elif timestep=='month':
      year=startyear
      month=startmonth-1
      for itime in range(length):
        month+=1
        timestr = str(year)+'-'+str(month)+'-01T00:00:00Z'
        dt_time = datetime.strptime(timestr, '%Y-%m-%dT%H:%M:%SZ')
        if unit=='hours since 1900-01-01 00:00:00':
          ini_timestr = str(1900)+'-'+str(1)+'-'+str(1)+'T00:00:00Z'
          ini_dt_time = datetime.strptime(ini_timestr, '%Y-%m-%dT%H:%M:%SZ')
          td = dt_time - ini_dt_time
          time_array.append(td.total_seconds()//3600)       
        if month==12:
          month=0
          year+=1
    return time_array    

def convert_101_1200_nc(nc_in_fn):
    print('Hello this this is convert_101_1200_nc')
    nc_dat = Dataset(nc_in_fn, 'r')
    varname = extract_varname.do(nc_dat)
    new_varname = varname.replace('101', '1200')
    basin_asc = ascraster.Asciigrid(ascii_file="/Users/pippo/Documents/SurfDrive/Research/Projects/DISCO_project/MyVersion/B_model_input/hydro/basin.map")
    dum_mon_nc = Dataset("/Users/pippo/Documents/SurfDrive/Research/Projects/DISCO_project/MyVersion/B_model_input/hydro/runoff_monthly.nc", "r")

    new_nc_in_fn = nc_in_fn.replace('101', '1200')
    new_nc_dat = Dataset(new_nc_in_fn, 'w')
    output_conversion.init_ncdata(os.getcwd(), new_nc_dat, new_varname, basin_asc)
    index = 0
    for itime in range(len(nc_dat['time'][0:100])):
      for i in range(12):
        grid=nc_dat[varname][itime,:,:]
        manip.add_grid_time(new_nc_dat, new_varname, grid, dum_mon_nc['time'][index])
        index +=1
    new_nc_dat.close()
    dum_mon_nc.close()
    nc_dat.close()
 
def convert_45_to_101_nc(nc_in_fn, interpolate='nearest'):
    print('Hello this this is convert_45_to_101_nc')
    nc_dat = Dataset(nc_in_fn, 'r')
    varname = extract_varname.do(nc_dat)
    time_array = make_time_array(101)
    basin_asc = ascraster.Asciigrid(ascii_file="/Users/pippo/Documents/SurfDrive/Research/Projects/DISCO_project/MyVersion/B_model_input/hydro/basin.map")
    new_nc_in_fn = nc_in_fn[:-3]+'_101.nc'
    new_nc_dat = Dataset(new_nc_in_fn, 'w')
    output_conversion.init_ncdata(os.getcwd(), new_nc_dat, varname, basin_asc,unit='', long_name='')
    print('entering the for loop')
    for itime in range(len(nc_dat['time'])):
        if itime < 14:
          for sub_itime in range(5):
              if interpolate==True:
                if not itime==(len(nc_dat['time'])-1):
                  try:
                    grid = nc_dat[varname][itime,:,:]*((5-(sub_itime))/5)+nc_dat[varname][itime+1,:,:]*((sub_itime)/5)
                  except(ZeroDivisionError):
                    grid = nc_dat[varname][itime,:,:]*((5-sub_itime)/5)+nc_dat[varname][itime+1,:,:]*0
                else:
                  grid = nc_dat[varname][itime,:,:]
              elif interpolate=='nearest':
                grid1 = nc_dat[varname][itime,:,:]
                grid2 = nc_dat[varname][itime+1,:,:]
                if sub_itime<1:
                  manip.add_grid_time(new_nc_dat, varname, grid1, time_array[sub_itime+(5*itime)])
                else:
                  manip.add_grid_time(new_nc_dat, varname, grid2, time_array[sub_itime+(5*itime)])
        else:
          grid = nc_dat[varname][itime,:,:]
          manip.add_grid_time(new_nc_dat, varname, grid, nc_dat['time'][itime])
    nc_dat.close()
    new_nc_dat.close()

def convert_27_to_101_nc(nc_in_fn, varname, interpolate=False):
    nc_dat = Dataset(nc_in_fn, 'r')
    time_array = make_time_array(101)
    basin_asc = ascraster.Asciigrid(ascii_file="/Users/pippo/Documents/SurfDrive/Research/Projects/DISCO_project/MyVersion/B_model_input/hydro/basin.map")
    new_nc_in_fn = nc_in_fn[:-3]+'_101.nc'
    new_nc_dat = Dataset(new_nc_in_fn, 'w')
    output_conversion.init_ncdata(os.getcwd(), new_nc_dat, varname+'_101', basin_asc,unit='[Mg C/km2/yr]', long_name='Aggregated net primary production (Mg C/km2/yr)')
    for itime in range(len(nc_dat['time'])):
        
        if itime < 7:
          if itime==0:
            for sub_itime in range(70):
              grid = nc_dat[varname][itime,:,:]
              manip.add_grid_time(new_nc_dat, varname+'_101', grid, time_array[sub_itime]) 
          else:
            for sub_itime in range(5):
              if interpolate:
                if not itime==(len(nc_dat['time'])-1):
                  try:
                    grid = nc_dat[varname][itime-1,:,:]*((5-(sub_itime))/5)+nc_dat[varname][itime,:,:]*((sub_itime)/5)
                  except(ZeroDivisionError):
                    grid = nc_dat[varname][itime-1,:,:]*((5-sub_itime)/5)+nc_dat[varname][itime,:,:]*0
                else:
                  grid = nc_dat[varname][itime,:,:] 
              else:
                grid = nc_dat[varname][itime,:,:] 
              manip.add_grid_time(new_nc_dat, varname+'_101', grid, time_array[sub_itime+(5*(itime-1))+70])
    manip.add_grid_time(new_nc_dat, varname+'_101', grid, time_array[100])
    nc_dat.close()
    new_nc_dat.close()

def do(args):
  if args[1]=='45_to_101':
    print('IF ELSE LOOP')
    convert_45_to_101_nc(args[2])
  elif args[1]=='101_to_1200':
    convert_101_1200_nc(args[2]) 
  elif args[1]=='27_to_101':
    convert_27_to_101_nc(args[2]) 

if __name__ == "__main__":
  do(sys.argv)
