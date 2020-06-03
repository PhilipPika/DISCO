# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import general python modules
import os
import sys
import traceback
import optparse
import numpy as np
from netCDF4 import Dataset
#from netCDF4 import MFDataset #read data from multiple NETCDF files at once
from netCDF4 import date2num,num2date,date2index
import time
import datetime as dt
import pickle

# From general module
import general_path
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
    #var_params = []                         #list of dictionaries
    
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
            #elif (type_option == '[SPECIE]'):
                #var_params.append(general_class.General())
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
                #elif (type_option == '[SPECIE]'):
                    #var_params[-1].add_item(option[0],option[1])
            except IndexError:
                raise MyError("File '%s' contains errors." % inifile,\
                              "Can not handle line: %s" % line)

    print('\n')
    # Close input file
    fp.close()

    # Check if all variables have a name
    # 2DO: other checks... src_params...
    #print '\n'
    #for ivar in range(len(var_params)):
        #if not('name' in var_params[ivar].get_attrib()):
            #raise MyError("File '%s' contains errors. No name was given to species nb %d." % (inifile,ivar+1))

    #return general_params, var_params
    #return general_params, src_params, var_params
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
    latlist = np.arange(yllcorner+nlat*cellsize-cellsize/2, yllcorner, -cellsize)
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

    
#def create_ncfile(ncfile, ncparams, varparams, nlat, nlon, iorder):
#    '''
#    Create NETCDF file for one species, with time, latitude and longitude dimensions
#    '''
#    rootgrp = Dataset(ncfile,'w',format=ncparams.get_val('format'))
#
#    # add try/except (error from general_class)
#    # Meta data of the file
#    rootgrp.description = ncparams.get_val('description') + ' for streams of order ' + str(iorder)
#    rootgrp.history = 'Created on ' + time.ctime(time.time())
#    rootgrp.source = ncparams.get_val('source')
#    
#    # Create dimensions: time is unlimited, others are fixed   
#    rootgrp.createDimension('time', None)
#    rootgrp.createDimension('lat', nlat)
#    rootgrp.createDimension('lon', nlon)
#
#    # Create variables 
#    date_time = rootgrp.createVariable('time','f4',('time',))
#    date_time.standard_name = 'time'
#    date_time.units = ncparams.get_val('time_units')
#    date_time.calendar = 'standard'
#    
#    lat = rootgrp.createVariable('lat','f4',('lat',))
#    lat.standard_name = 'latitude'
#    lat.units = 'degrees north'
#    
#    lon = rootgrp.createVariable('lon','f4',('lon',))
#    lon.standard_name = 'longitude'
#    lon.units = 'degrees east'
#
#    add_variable(rootgrp, ncparams, varparams)
#    
#    rootgrp.sync()
#    print rootgrp
#    rootgrp.close()


#def create_ncfiles(spec_orders, general_params, var_params, nlat, nlon):
#    '''
#    Initializes NETCDF files based on description from inifile
#    Creates files for every variable and each order it is distributed to
#    spec_orders is a dictionary of lists of orders to which each species is distributed
#    Returns a dictionary with species nb as key of [order, filename]
#    '''
#    filenames = {}
#
#    for ivar in range(len(var_params)):
#        filenames[ivar] = []
#        
#    # Make new output directory when it does not exist
#    create_dir(general_params.get_val('path'))
#
#    for ivar in range(len(var_params)):
#        for iorder in spec_orders[ivar]:
#            filename = os.path.join(general_params.get_val('path'), var_params[ivar].get_val('name')+'_'+str(iorder)+'.nc')
#            create_ncfile(filename, general_params, var_params[ivar], nlat, nlon, order=iorder)
#            filenames[ivar].append([iorder, filename])
#
#    return filenames


#def add_grid_time(ncfile, varname, grid, t0):
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


def find_within_range(times, starttime, endtime, ncfile):
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
    for key in rootgrp.variables.keys():
        if (key.upper() == varname.upper()):
            return key
    else:
        print('***** ERROR *****')
        print('Specified variable ' + varname + ' is not found.')
        print('Names available: ')
        for key in rootgrp.variables.keys():
            print(key)
        raise MyError("File '%s' does not exist" % ncfile)


###################
### Conversions ###
###################
    
def create_ncmask_from_ascii_file(asciifile):
    '''
    Returns a mask (numpy array of booleans) corresponding to ASCII file
    Grid cells with no data values in ASCII files are masked
    And a dictionnary containing properties of the ASCII grid (nb of rows, cols, xllcorner, yllcorner and cellsize)
    '''

    # Open ASCII file
    ascii_map = ascraster.Asciigrid(ascii_file=asciifile)

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
    for icell in xrange(ascii_map.length):
        #print 'data ' + str(icell) + ': ' + str(map.get_data(icell))
        ilat,ilon = ascii_map.get_row_col_from_index(icell)

        #print 'ASCII map value for cell ' + str(icell) + '  ' + str(ascii_map.get_data(icell))
        if (ascii_map.get_data(icell) == None):
            netcdf_mask[ilat, ilon] = True

        #print 'netcdf_mask  ' + str(netcdf_mask)

    return netcdf_mask, ascii_map_prop


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

def calc_row_col_from_index(ind, ncols):
    '''
    Returns irow and icol of an index in an unmasked grid
    '''
    irow = int(int(ind)/int(ncols))
    icol = int(ind) - irow * ncols
    return irow,icol


def calc_index_from_row_col(irow, icol, nrows, ncols):
    '''
    Get index of combination of row and col numbers
    '''
    ind = irow * int(ncols) + icol
    grid_length = nrows * ncols
    #print 'ind  ' + str(ind)
    if (ind >= grid_length) or (ind < 0): 
        raise MyError("calc_index_from_row_col: Row,Col position (%s,%s) falls outside bounds." % (irow,icol))
    return ind

def calc_index_from_coordin(x, y, nrows, ncols, xllcorner, yllcorner, cellsize):
    '''
    Get index of x- and y-coordinate
    '''
    irow = nrows - int((float(y)-float(yllcorner))/float(cellsize)) - 1
    #print 'irow  ' + str(irow)
    icol = int((float(x)-float(xllcorner))/float(cellsize))
    #print 'icol  ' + str(icol)
    
    if (icol >= ncols or icol < 0):
        raise MyError("calc_index_from_coordin: Y position (%s,%s) falls outside bounds." % (x,y))
    if (irow >= nrows or irow < 0):
        raise MyError("calc_index_from_coordin: X position (%s,%s) falls outside bounds." % (x,y))
    
    return calc_index_from_row_col(irow, icol, nrows, ncols)
    

def calc_row_col_from_coordin(x, y, nrows, ncols, xllcorner, yllcorner, cellsize):
    '''
    Returns irow and icol from coordinates
    '''
    ind = calc_index_from_coordin(x, y, nrows, ncols, xllcorner, yllcorner, cellsize)
    return calc_row_col_from_index(ind, ncols)

    
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
    
    #ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31]
    #dtdates = num2date(numdates, units)
    #years = []

    # num2date returns a single datetime object if only one time is to be converted
    # and a numpy array of datetime objects if more than one need conversion
    #if (type(dtdates) == dt.datetime):
    #    dtdates = [dtdates]
    
    #for idate in range(len(dtdates)):
    #    year = float(dtdates[idate].year)
    #    ndays = 365
    #    if (year%400==0 or (year%4==0 and year%100!=0)):
    #        ndays += 1
    #    days = sum(ndays_mon[0:dtdates[idate].month-1])
    #    days += dtdates[idate].day-1 #We consider the value is always given at hour 00:00:00 - 2 change
    #    year += float(days) / ndays
    #    years.append(year)

    #if (len(years) == 1):
    #    return years[0]
    #else:
    #    return years
    

############
### Test ###
############
if __name__ == "__main__":
    
    try:
        if (len(sys.argv) < 2):
            raise MyError("Not enough arguments on command line.\n" +\
                          "Should be: python manip_netcdf.py name_inifile")
        create_ncfiles(sys.argv[1])
    #except MyError, val:
        #val.write()
    #except (optparse.OptionValueError,
    #        Exception), val:
    #    print(str(val))
    #    print(sys.exc_type)
    #    print(traceback.print_exc())
    except:
        print('***** ERROR ******')
        print('manip.py failed.')
        print(sys.exc_type)
        print(traceback.print_exc())
