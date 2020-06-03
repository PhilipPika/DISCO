# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/manip_netcdf.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import general python modules
import os
import sys
import traceback
import optparse
import numpy as np
from netCDF4 import Dataset
from netCDF4 import date2num,num2date,date2index
import time
import datetime as dt
import pickle as pickle
import math

# From general module
import ascraster
import general_class
from error import *

# Raise division by zero errors when using numpy
np.seterr(divide='raise')


#############################################################
### Create NETCDF files with parameters read in a inifile ###
#############################################################
def read_parameters_ncinifile(inifile):
    '''
    Reads general and variable parameters in an inifile
    Stores it and returns in a dictionary form
    '''
    general_params = general_class.General() #dictionary
    src_params = []                          #list of dictionaries
    
    check_file(inifile)
  
    # Open file with all parameters
    fp = open(inifile,'r')

    # Determine all the information per line
    type_option = None
    for line in fp.readlines():
        option = [elem.strip() for elem in line.split('#')[0].split('=')]
        if (len(option[0]) < 1):
            # Empty line
            continue
        
        # Check whether the first element is start of a new class 
        # (first character is a '[' )
        if option[0].startswith('['):
            # Start of a new class
            type_option = option[0].strip().upper()
            print('\n' + type_option)
            if (type_option == '[SOURCE]'):
                src_params.append(general_class.General())
            elif (type_option != '[GENERAL]'):
                raise MyError("File '%s' contains errors. Class %s is not known." % (inifile,option[0]))

        else:
            # Here information from a class is given
            print(option[0] + ' = ' + option[1])
            try:
                if (type_option == '[GENERAL]'):
                    general_params.add_item(option[0],option[1])
                elif (type_option == '[SOURCE]'):
                    src_params[-1].add_item(option[0],option[1])
            except IndexError:
                raise MyError("File '%s' contains errors." % inifile,\
                              "Can not handle line: %s" % line)
    # Close input file
    fp.close()

    return general_params, src_params


def add_variable(rootgrp, ncparams, varparams):
    '''
    Add a new variable in NETCDF dataset
    Fill in variable metadata
    '''
    varname = varparams.get_val('name')    
    var = rootgrp.createVariable(varname,'f4',('time','lat','lon',) ,fill_value=ncparams.get_val('fill'),zlib=ncparams.get_val('zlib'))
    var.standard_name = varname
    if ('units' in varparams.get_attrib()):
        var.units = varparams.get_val('units')
    if ('description' in varparams.get_attrib()):
        var.description = varparams.get_val('description')

    
def initialize_lat(nlat, yllcorner, cellsize):
    '''
    Calculates latitudes from yllcorner value and cellsize
    '''
    latlist = np.arange(yllcorner+nlat*cellsize-cellsize/2., yllcorner, -cellsize)
    return latlist
    
    
def initialize_lon(nlon, xllcorner, cellsize):
    '''
    Calculates longitudes from xllcorner value and cellsize
    '''
    lonlist = np.arange(xllcorner+cellsize/2., xllcorner+nlon*cellsize, cellsize)
    return lonlist
    
    
def create_time_lat_lon(rootgrp, lat, lon, time_units, time_calendar='standard'):
    '''
    Create time, latitude and longitude dimensions and variables in new NETCDF datset
    '''

    nlat = len(lat)
    nlon = len(lon)
    
    # Create dimensions: time is unlimited, others are fixed   
    rootgrp.createDimension('time', None)
    rootgrp.createDimension('lat', nlat)
    rootgrp.createDimension('lon', nlon)

    # Create variables 
    date_time = rootgrp.createVariable('time','f4',('time',))
    date_time.standard_name = 'time'
    date_time.units = time_units
    date_time.calendar = time_calendar
    
    latitude = rootgrp.createVariable('lat','f4',('lat',))
    latitude.standard_name = 'latitude'
    latitude.units = 'degrees north'
    latitude[:] = lat
    
    longitude = rootgrp.createVariable('lon','f4',('lon',))
    longitude.standard_name = 'longitude'
    longitude.units = 'degrees east'
    longitude[:] = lon

def add_grid_time(ncdata, varname, grid, t0):
    '''
    Add a grid of values for of values for one variable at time t0 in NETCDF file
    '''
    #ncdata = Dataset(ncfile, 'a')
    #print 't0  ' + str(t0)
    nt = len(ncdata.variables['time'][:])
    ncdata.variables['time'][nt] = t0
    ncdata.variables[varname][nt,:,:] = grid
    #ncdata.close() #inefficient to open/close to add each grid?


################################################
### Read and extract data from a NETCDF file ###
################################################
def get_netcdf_year(netcdf_filename,varname,year,mask=None,nodata_value=-9999.):
    '''
    Returns the whole dataset at a given time.
    If prescribed (ncols,nrows,xll,yll and cellsize in params) then that format
    is returned. Also if mask is specified then masked grid is returned.
    If params is None, then netcdf as is is returned for that year.
    Nodata_value is used as default nodata value of the output grid.
    idim is the number of the fourth dimension. Only when the number of dimensions is four,
    you must specify the last dimension. Normal {time,lat,lon,ndim} as dimensions for this.
    '''
    # We have to open and close the filename
    if (not os.path.isfile(netcdf_filename)):
        raise MyError("File: " + netcdf_filename + " does not exist.")
    ncdata = Dataset(netcdf_filename, 'r')

    # Determine the years that are given for this parameter
    nctimes = ncdata.variables['time']
    years = convert_numdate2year(nctimes[:], nctimes.units)
    
    # Find begin and end point for this time period
    istart,iend = find_within_range(years,year,year)

    # Check whether one of the boundaries is equal to the year specified.
    # Range of 0.001 is used, because this is smaller than 0.5 day.
    if (math.fabs(years[istart] - year) < 0.001):
        ifound = istart
    elif (math.fabs(years[iend-1] - year) < 0.001):
        ifound = iend-1
    else:
        raise MyError("File: " + netcdf_filename + " does not have data for: " + str(year) +".",\
                      "Years found: " + ";".join(map(str,years)))

    # Assign the two dimension data for this parameter.
    try:
        vardata = ncdata.variables[varname][ifound,:,:]
    except KeyError:
        raise MyError("File: " + netcdf_filename + " does not contain parameter: " + varname,\
                      "Parameters available: " + ";".join(map(str,list(ncdata.variables.keys())))+".")

    # Create a Asciigrid of this data
    # Get latitude
    try:
        lat = ncdata.variables['latitude'][:]
    except KeyError:
        try:
            lat = ncdata.variables['lat'][:]
        except KeyError:
            raise MyError("File: " + netcdf_filename + " does not contain parameter: lat or latitude.",\
                          "Parameters available: " + ";".join(map(str,list(ncdata.variables.keys())))+".")

    # Set number of rows        
    nrows = len(lat)

    # Get longitude
    try:
        lon = ncdata.variables['longitude'][:]
    except KeyError:
        try:
            lon = ncdata.variables['lon'][:]
        except KeyError:
            raise MyError("File: " + netcdf_filename + " does not contain parameter: lon or longitude.",\
                          "Parameters available: " + ";".join(map(str,list(ncdata.variables.keys())))+".")
    # Set number of columns
    ncols = len(lon)

    # Get stepsize of the grid.
    stepsize_col = math.fabs((lon[-1]-lon[0])/float(ncols-1))
    stepsize_row = math.fabs((lat[-1]-lat[0])/float(nrows-1))

    if (math.fabs(stepsize_col - stepsize_row) > 0.0001):
        raise MyError("File: " + netcdf_filename + " has not the same cellsize in lat - lon direction.")

    # Set cellsize and xllcorner and yllcorner
    cellsize = stepsize_col
    xllcorner = min(lon[-1],lon[0]) - 0.5 * cellsize
    yllcorner = min(lat[-1],lat[0]) - 0.5 * cellsize

    # Find out whether the data must be flipped.
    # We want to start left-upper corner and then go rowwise down.
    if (lat[0] < lat[-1]):
        # We have to flip (bottom to top and top to bottom)
        vardata_flip = np.flipud(vardata)
        vardata = vardata_flip
        print("Flipud is performed. Must be checked!")
    if (lon[-1] < lon[0]):
        # We have to flip (left-right)
        vardata_flip = np.fliplr(vardata)
        vardata = vardata_flip
        print("Fliplr is performed. Must be checked!")

    # Reshape the three dimension grid into two dimensional list.
    try:
        vardata_list = list(vardata.filled([nodata_value]).reshape((nrows*ncols,-1)))
    except AttributeError:
        # Array is not masked.
        vardata_list = list(vardata.reshape((nrows*ncols,-1)))

    # Apply the mask when it is given
    if (mask != None):
        vardata_list = [vardata_list[int(x)] for x in mask]

    # Close file when we opened it here.
    ncdata.close() 

    # Return years
    return vardata_list


def get_netcdf_year_grid(netcdf_filename,varname,year,params=None,mask=None,nodata_value=-9999.,idim=0):
    '''
    Returns the whole grid at a given time.
    If prescribed (ncols,nrows,xll,yll and cellsize in params) then that format
    is returned. Also if mask is specified then masked grid is returned.
    If params is None, then netcdf as is is returned for that year.
    Nodata_value is used as default nodata value of the output grid.
    idim is the number of the fourth dimension. Only when the number of dimensions is four,
    you must specify the last dimension. Normal {time,lat,lon,ndim} as dimensions for this.
    '''
    # We have to open and close the filename
    if (not os.path.isfile(netcdf_filename)):
        raise MyError("File: " + netcdf_filename + " does not exist.")
    ncdata = Dataset(netcdf_filename, 'r')

    # Determine the years that are given for this parameter
    nctimes = ncdata.variables['time']
    years = convert_numdate2year(nctimes[:], nctimes.units)
    
    # Find begin and end point for this time period
    istart,iend = find_within_range(years,year,year)

    # Check whether one of the boundaries is equal to the year specified.
    # Range of 0.001 is used, because this is smaller than 0.5 day.
    if (math.fabs(years[istart] - year) < 0.001):
        ifound = istart
    elif (math.fabs(years[iend-1] - year) < 0.001):
        ifound = iend-1
    else:
        raise MyError("File: " + netcdf_filename + " does not have data for: " + str(year) +".",\
                      "Years found: " + ";".join(map(str,years)))

    # Assign the two dimension data for this parameter.
    try:
        if (len(ncdata.variables[varname].dimensions) == 3):
            vardata = ncdata.variables[varname][ifound,:,:]

        elif (len(ncdata.variables[varname].dimensions) == 4):
            try:
                vardata = ncdata.variables[varname][ifound,:,:,idim]
            except IndexError:
                raise MyError("File: " + netcdf_filename + " contains parameter: " + varname,\
                              "Dimension given (given with idim): " + str(idim) + " is not available.",\
                              "Available dimension must be between 0 and " + \
                              str(len(ncdata.variables[varname][ifound,0,0,:])-1))
 
        else:
            raise MyError("File: " + netcdf_filename + " contains parameter: " + varname + \
                          " with unknown dimensions.")
    except KeyError:
        raise MyError("File: " + netcdf_filename + " does not contain parameter: " + varname,\
                      "Parameters available: " + ";".join(map(str,list(ncdata.variables.keys())))+".")

    # Create a Asciigrid of this data
    # Get latitude
    try:
        lat = ncdata.variables['latitude'][:]
    except KeyError:
        try:
            lat = ncdata.variables['lat'][:]
        except KeyError:
            raise MyError("File: " + netcdf_filename + " does not contain parameter: lat or latitude.",\
                          "Parameters available: " + ";".join(map(str,list(ncdata.variables.keys())))+".")

    # Set number of rows        
    nrows = len(lat)

    # Get longitude
    try:
        lon = ncdata.variables['longitude'][:]
    except KeyError:
        try:
            lon = ncdata.variables['lon'][:]
        except KeyError:
            raise MyError("File: " + netcdf_filename + " does not contain parameter: lon or longitude.",\
                          "Parameters available: " + ";".join(map(str,list(ncdata.variables.keys())))+".")
    # Set number of columns
    ncols = len(lon)

    # Get stepsize of the grid.
    stepsize_col = math.fabs((lon[-1]-lon[0])/float(ncols-1))
    stepsize_row = math.fabs((lat[-1]-lat[0])/float(nrows-1))

    if (math.fabs(stepsize_col - stepsize_row) > 0.0001):
        raise MyError("File: " + netcdf_filename + " has not the same cellsize in lat - lon direction.")

    # Set cellsize and xllcorner and yllcorner
    cellsize = stepsize_col
    xllcorner = min(lon[-1],lon[0]) - 0.5 * cellsize
    yllcorner = min(lat[-1],lat[0]) - 0.5 * cellsize

    # Find out whether the data must be flipped.
    # We want to start left-upper corner and then go rowwise down.
    if (lat[0] < lat[-1]):
        # We have to flip (bottom to top and top to bottom)
        vardata_flip = np.flipud(vardata)
        vardata = vardata_flip
        print("Flipud is performed.")
    if (lon[-1] < lon[0]):
        # We have to flip (left-right)
        vardata_flip = np.fliplr(vardata)
        vardata = vardata_flip
        print("Fliplr is performed.")

    # Set nodata value to the specified value
    #vardata1 = vardata.filled([nodata_value])
    #vardata = vardata1
    # Reshape the two dimension grid into one dimensional list.
    try:
        vardata_list = list(np.reshape(vardata.filled([nodata_value]), nrows*ncols))
    except AttributeError:
        # Array is not masked.
        vardata_list = list(np.reshape(vardata, nrows*ncols))

    # Create output grid.
    if (params == None):
        outgrid = ascraster.Asciigrid(nrows=nrows,ncols=ncols,xllcorner=xllcorner,\
                                      yllcorner=yllcorner,cellsize=cellsize,\
                                      nodata_value = nodata_value)

    else:
        # Dependent on cellsize, choose the right parameters setting.
        if (math.fabs(cellsize - 0.5) < 0.0001):
            # Half degree settings.
            try:
                outgrid = ascraster.Asciigrid(nrows=params.nrows,ncols=params.ncols,xllcorner=params.xllcorner,\
                                          yllcorner=params.yllcorner,cellsize=params.cellsize,\
                                          nodata_value = nodata_value)
            except AttributeError:
                raise MyError("Set the attributes nrows, ncols, xllcorner, yllcorner and cellsize",\
                          "in the class params for the function manip_netcdf.get_netcdf_year_grid.")

        elif(math.fabs(cellsize - 0.08333) < 0.0001):
            # Five minutes settings.
            try:
                outgrid = ascraster.Asciigrid(nrows=params.nrows_5min,ncols=params.ncols_5min,xllcorner=params.xllcorner,\
                                          yllcorner=params.yllcorner,cellsize=params.cellsize_5min,\
                                          nodata_value = nodata_value)
            except AttributeError:
                raise MyError("Set the attributes nrows_5min, ncols_5min, xllcorner, yllcorner and cellsize_5min",\
                          "in the class params for the function manip_netcdf.get_netcdf_year_grid.")
        else:
            raise MyError("File: " + netcdf_filename + " does not have the expected spatial resolution.",\
                          "Read this file with params=None in function manip_netcdf.get_netcdf_year_grid.")

    # Put data into grid.
    outgrid.add_values(vardata_list)

    # Apply the mask when it is given
    if (mask != None):
        outgrid.apply_mask(mask)

    # Close file when we opened it here.
    ncdata.close() 

    # Return years
    return outgrid



def get_value_at_date(rootgrp, varname, t0, ilat, ilon, no_data_value=None):
    '''
    Returns the value of a variable contained in a NETCDF file at a given time in cell (ilat, ilon)
    If the time is not available in NETCDF file, interpolation is performed (with extrapolation of 1st/last available value)
    '''
    # Check whether variable name is available
    key = check_variable_in_ncfile(rootgrp, varname, ncfile)

    val = None
    #print 'Getting value at (ilat, ilon) = (' + str(ilat) + ', ' + str(ilon) + ').'
    # Check if values are missing in the cell
    if (type(rootgrp.variables[key][:,ilat,ilon]) is np.ma.masked_array):
        print('***** WARNING *****')
        print('No value for ' + key + ' at (ilat, ilon) = (' +\
            str(ilat) + ', ' + str(ilon) + ').')
        print('Value set to ' + str(no_data_value))
        val = no_data_value
        
    # If not, extraction of the value t time t0
    else:
        #2DO: add check for None values
        val = np.interp(t0, rootgrp.variables['time'][:],\
                        rootgrp.variables[key][:,ilat,ilon])

    return val


def find_within_range_netcdf(times, starttime, endtime, ncfile):
    '''
    Returns the indexes of the times in ncfile that correspond to the beginning and end of the simulation period
    '''
    if (starttime < times[0]):
        raise MyError("Start time is inferior to available times in '%s'" % ncfile)
    if (endtime < times[-1]):
        raise MyError("End time is superior to available times in '%s'" % ncfile)
    istart = date2index(starttime, times, select='before')
    iend = date2index(endtime, times, select='after')
    return istart,iend

def find_within_range(list1,begin,end):
    '''
    Find all the values which are needed to cover the range begin to end.
    This means that the last element before begin is taken into acount and the the first after end is included.
    list1 must be sorted!
    '''
    ibegin = 0
    iend = len(list1) 
    for item in range(1,len(list1)):
        if (float(list1[item]) < begin):
            ibegin = item
    for item in range(len(list1)-1,-1,-1):
        if (float(list1[item]) > end):
            iend = item
    return ibegin,iend

def get_times(starttime=None,endtime=None,ncdata=None,netcdf_filename=None):
    '''
    Returns the times that are given for this parameter in a list.
    Or ncdata (Dataset) is given or the filename must be given.
    '''
    
    if (ncdata == None):
        if (netcdf_filename == None):
            raise MyError("Wrong call to function manip_netcdf.get_times.")
        else:
            # We have to open and close the filename
            if (not os.path.isfile(netcdf_filename)):
                raise MyError("File: " + netcdf_filename + " does not exist.")
            ncdata = Dataset(netcdf_filename, 'r')

    # Determine the years that are given for this parameter
    nctimes = ncdata.variables['time']
    years = convert_numdate2year(nctimes[:], nctimes.units)
    
    # Find begin and end point for this time period
    if (starttime != None and endtime != None):
        istart,iend = find_within_range(years,starttime,endtime)
        years = years[istart:iend]

    # Close file when we opened it here.
    if (ncdata == None):
        ncdata.close()     

    # Return years
    return years

##########################################
### Convert ASCII files to NETCDF file ###
##########################################


##############
### Checks ###
##############
def create_dir(dirname):
    '''
    Creates  directory if it doesn't exist yet
    '''
    if not (os.path.isdir(dirname)):
        os.mkdir(dirname)

def check_file(filename):
    '''
    Checks if file exists
    '''
    if not os.path.exists(filename):
        raise MyError("File '%s' does not exist" % filename)

def check_variable_in_ncfile(rootgrp, varname, ncfile):
    '''
    Checks is a variable exists in NETCDF file
    '''
    for key in list(rootgrp.variables.keys()):
        if (key.upper() == varname.upper()):
            return key
    else:
        print('***** ERROR *****')
        print('Specified variable ' + varname + ' is not found.')
        print('Names available: ')
        for key in list(rootgrp.variables.keys()):
            print(key)      
        raise MyError("File '%s' does not exist" % ncfile)


###################
### Conversions ###
###################
    
def create_ncmask_from_ascii_file(asciifile, numtype=float):
    '''
    Returns a mask (numpy array of booleans) corresponding to ASCII file
    Grid cells with no data values in ASCII files are masked
    And a dictionnary containing properties of the ASCII grid (nb of rows, cols, xllcorner, yllcorner and cellsize)
    '''

    # Open ASCII file
    ascii_map = ascraster.Asciigrid(ascii_file=asciifile, numtype=numtype)

    # Create empty mask (nothing is masked yet)
    netcdf_mask = np.zeros((ascii_map.nrows, ascii_map.ncols), dtype=bool)

    # Create dictionnary containing map's properties
    ascii_map_prop = {}
    ascii_map_prop['ncols'] = ascii_map.ncols
    ascii_map_prop['nrows'] = ascii_map.nrows
    ascii_map_prop['xllcorner'] = ascii_map.xllcorner
    ascii_map_prop['yllcorner'] = ascii_map.yllcorner
    ascii_map_prop['cellsize'] = ascii_map.cellsize

    # For each cell, mask if no data in the ASCII grid
    #print 'Length of ASCII map  ' + str(ascii_map.length)
    #print 'ASCII map nodata value  ' + str(ascii_map.nodata_value)
    for icell in range(ascii_map.length):
        #print 'data ' + str(icell) + ': ' + str(map.get_data(icell))
        ilat,ilon = ascii_map.get_row_col_from_index(icell)

        #print 'ASCII map value for cell ' + str(icell) + '  ' + str(ascii_map.get_data(icell))
        if (ascii_map.get_data(icell) == None):
            netcdf_mask[ilat, ilon] = True

        #print 'netcdf_mask  ' + str(netcdf_mask)

    return netcdf_mask, ascii_map_prop


def ncgrid_from_asciifile(asciifile, nlat, nlon, netcdf_mask, fill_value=None):
    '''
    Creates an masked numpy array grid (to be included in NETCDF file)
    from an ASCII file
    '''

    # Open ASCII grid
    ascii_grid = ascraster.Asciigrid(ascii_file=asciifile, numtype=float)

    # Create empty masked numpy array
    netcdf_grid = np.ma.array(np.zeros((nlat, nlon)), mask=netcdf_mask, fill_value=fill_value)

    # Fill numpy array with values from the ASCII grid
    #for item in range(ascii_grid.length):
        #ilat,ilon = ascii_grid.get_row_col_from_index(item)
    for ilat in range(nlat):
        for ilon in range(nlon):
            item = ascii_grid.get_index_from_row_col(ilat, ilon)
            if not(netcdf_mask[ilat,ilon]): #TEST
                if (ascii_grid.get_data(item) == None):
                    print('Oops! None value replaced by zero')
                    netcdf_grid[ilat,ilon] = 0.
                else:
                    netcdf_grid[ilat,ilon] = ascii_grid.get_data(item)
    
    return netcdf_grid


def calc_lat_lon_from_mask(mask, nlat, nlon):
    '''
    Calculates the indexes of latitude and longitude dimensions from an ASCII mask
    '''
    # Use in rive_riverbasin when creating the rivmask?
    lat = []
    lon = []

    for item in range(len(mask)):
        icell = mask[item]
        ilat = int(icell/nlon)
        ilon = icell - ilat*nlon
        lat.append(ilat)
        lon.append(ilon)
        
    return lat,lon

def calc_row_col_from_index(ind,ncols):
    '''
    Returns irow and icol of an index
    '''
    if (type(ncols) != int or type(ind) != int):
        raise MyError("Function calc_row_col_from_index must have integers as argument")
    # Note, here a division between integers is needed. 
    irow = ind//ncols
    icol = ind - irow * ncols
    return irow,icol


def calc_index_from_row_col(irow, icol, nrows, ncols):
    '''
    Get index of combination of row and col numbers
    '''
    if (type(ncols) != int or type(nrows) != int or type(icol) != int or type(irow) != int):
        raise MyError("Function calc_index_from_row_col must have integers as argument")

    return irow * ncols + icol

def calc_index_from_coordin(x, y, nrows, ncols, xllcorner, yllcorner, cellsize):
    '''
    Get index of x- and y-coordinate
    '''
    irow, icol = calc_row_col_from_coordin(x, y, nrows, ncols, xllcorner, yllcorner, cellsize)
    
    return calc_index_from_row_col(irow, icol, nrows, ncols)
    

def calc_row_col_from_coordin(x, y, nrows, ncols, xllcorner, yllcorner, cellsize):
    '''
    Returns irow and icol from coordinates
    '''
    irow = nrows - int((float(y) - float(yllcorner))/float(cellsize)) - 1
    icol = int((float(x) - float(xllcorner))/float(cellsize))
    
    if (icol >= ncols or icol < 0):
        raise MyError("calc_index_from_coordin: Y position (%s,%s) falls outside bounds." % (x,y))
    if (irow >= nrows or irow < 0):
        raise MyError("calc_index_from_coordin: X position (%s,%s) falls outside bounds." % (x,y))
    
    return irow, icol

    
def convert_year2numdate(year0, units):
    '''
    Converts a date in years (decimal) to number used in NETCDF file
    units should be in format <unit> since yyyy-mm-dd hh:mm:ss
    <unit> can be days, hours, ...
    '''
    date_tmp = convert_year2datetime(year0)
    return date2num(date_tmp, units)


def convert_year2datetime(year0):
    '''
    Converts a date in years (decimal) to a datetime object
    '''
    year = int(year0)
    ndays = 365
    if (year%400==0 or (year%4==0 and year%100!=0)):
        ndays += 1
    rest = (year0-year) * ndays
    units_tmp = 'days since ' + str(year) + '-01-01 00:00:00'
    return num2date(rest, units_tmp)

def convert_datetime2year(dtdates):
    '''
    Converts a list of datetime objects to decimal years
    '''
    ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31]
    years = []

    # num2date returns a single datetime object if only one time is to be converted
    # and a numpy array of datetime objects if more than one need conversion
    if (type(dtdates) == dt.datetime):
        dtdates = [dtdates]
    
    for idate in range(len(dtdates)):
        year = float(dtdates[idate].year)
        ndays = 365
        if (year%400==0 or (year%4==0 and year%100!=0)):
            ndays += 1
        days = sum(ndays_mon[0:dtdates[idate].month-1])
        days += dtdates[idate].day-1 #We consider the value is always given at hour 00:00:00 - 2 change
        year += float(days) / ndays
        years.append(year)

    if (len(years) == 1):
        return years[0]
    else:
        return years

def convert_numdate2year(numdates, units):
    '''
    Converts a list of times in numbers (units defined in NETCDF file) to decimal years
    '''
    dtdates = num2date(numdates, units)
    years = convert_datetime2year(dtdates)
    return years
    

############
### Test ###
############
if __name__ == "__main__":
    
    try:
        #create_ncfiles(sys.argv[1])

        netcdf_filename="/home/arthur/SSP2/netcdf/exchange/GFRAC.nc"
        mask = list(range(1043503,1063504))
        year = 2000. 
        dataset = get_netcdf_year(netcdf_filename,"GFRAC",year,mask=mask,nodata_value=-9999.)
        print(len(dataset),len(dataset[0]))
        for icell in range(len(dataset)):
            if (sum(dataset[icell]) > 0):
                print(icell,year,dataset[icell],sum(dataset[icell]))
        year = 2005. 
        dataset = get_netcdf_year(netcdf_filename,"GFRAC",year,mask=mask,nodata_value=-9999.)
        print(len(dataset),len(dataset[0]))
        for icell in range(len(dataset)):
            if (sum(dataset[icell]) > 0):
                print(icell,year,dataset[icell],sum(dataset[icell]))
        year = 2010. 
        dataset = get_netcdf_year(netcdf_filename,"GFRAC",year,mask=mask,nodata_value=-9999.)
        print(len(dataset),len(dataset[0]))
        for icell in range(len(dataset)):
            if (sum(dataset[icell]) > 0):
                print(icell,year,dataset[icell],sum(dataset[icell]))

        #raise MyError("KLAAR")
        print(get_times(starttime=2000,endtime=2045,netcdf_filename=netcdf_filename))
        print("***")
        print(get_times(starttime=2000.,endtime=2000.,netcdf_filename=netcdf_filename))
        print("***")
        print(get_times(netcdf_filename=netcdf_filename))
        print("***")
        ncdata = Dataset(netcdf_filename,'r')
        print(get_times(ncdata=ncdata))
        ncdata.close()
        import cmd_options_general
        params = cmd_options_general.Input()
        params.nrows_5min = 2160
        params.ncols_5min = 4320
        params.xllcorner = -180
        params.yllcorner = -90
        params.cellsize_5min = 0.083333
        try:
            netcdf_filename="/home/arthur/SSP2/netcdf/exchange/GFRAC.nc"
            outgrid = get_netcdf_year_grid(netcdf_filename,"GFRAC",2001.,nodata_value=-1)
            outgrid.write_ascii_file("piet.asc")
        except MyError as val:
            val.write()
        try:
            netcdf_filename="/home/arthur/SSP2/netcdf/exchange/GFRAC.nc"
            outgrid = get_netcdf_year_grid(netcdf_filename,"gfrac",2000.,nodata_value=-1)
            outgrid.write_ascii_file("piet.asc")
        except MyError as val:
            val.write()
        try:
            netcdf_filename="/home/arthur/SSP2/netcdf/exchange/GFRAC.nc"
            outgrid = get_netcdf_year_grid(netcdf_filename,"GFRAC",2000.,idim=20,nodata_value=-1)
            outgrid.write_ascii_file("piet.asc")
        except MyError as val:
            val.write()
        try:
            netcdf_filename="/home/arthur/SSP2/netcdf/exchange/GFR.nc"
            outgrid = get_netcdf_year_grid(netcdf_filename,"GFRAC",2000.,nodata_value=-1)
            outgrid.write_ascii_file("piet.asc")
        except MyError as val:
            val.write()
        #raise MyError("klaar")
        netcdf_filename="/home/arthur/SSP2/netcdf/exchange/GNDEP.nc"
        outgrid = get_netcdf_year_grid(netcdf_filename,"GNDEP",2000.,params,idim=1,nodata_value = -2)
        outgrid.write_ascii_file("piet_dep.asc")
        netcdf_filename="/home/arthur/SSP2/netcdf/exchange/GNDEP.nc"
        outgrid = get_netcdf_year_grid(netcdf_filename,"GNDEP",2000.)
        outgrid.write_ascii_file("piet_dep1.asc")

    except MyError as val:
        val.write()
    except (optparse.OptionValueError,
            Exception) as val:
        print(str(val))
        print(sys.exc_info()[0])
        print(traceback.print_exc())
    except:
        print('***** ERROR ******')
        print('manip_netcdf.py failed.')
        print(sys.exc_info()[0])
        print(traceback.print_exc())
