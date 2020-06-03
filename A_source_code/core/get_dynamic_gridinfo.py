# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import Python modules
import numpy as np
from netCDF4 import Dataset, MFDataset
import os
import sys
import traceback
import optparse
import pickle
import zipfile
import copy

# Import own general modules
import manip_netcdf
import pointer_class
import general_class
import directory
from error import *

# Import local modules.
import interpolate_list
import general_func


# For test
class Input_test():
    def __init__(self):
       pass


def calculate_yearly_load(load_cell, year_start, year_end):
    '''
    Calculates a yearly time series from year_start to year_end (should be integers)
    '''
    #year_start = int(load_cell[0][0])
    #year_end = int(load_cell[-1][0])
    load_yearly = []

    # The value attributed to the first year is the first available value
    #load_yearly.append(load_cell[0])
    
    for year in range(year_start, year_end+1):
        val = interpolate_list.calculate(year, load_cell, extrapol=1)
        load_yearly.append([year, val[1]])
        
    return load_yearly


def calculate_disaggregated_load(params, load_yearly, ipointer, temp_distrib):
    '''
    Applies a temporal pattern to a yearly time series
    The returned time series has the same time steps as the provided pattern
    '''
    year_start = int(load_yearly[0][0])
    year_end = int(load_yearly[-1][0])+1 #calculations stop at the end of year_end
    load_dis = []

    pattern_file = '' #file containing the temporal pattern variable
    varname = '' #name of the temporal pattern variable
    square = False

    # Open the pattern file and define the variable name, depending on the provided temp_distrib label
    if (temp_distrib == 'uniform'):
        pattern_file = os.path.join(params.water_inputdir,params.uniform)
        varname = params.uniform_varname
    elif (temp_distrib == 'discharge'):
        pattern_file = os.path.join(params.water_inputdir,params.discharge)
        varname = params.discharge_varname
    elif (temp_distrib == 'flooding_volume'):
        pattern_file = os.path.join(params.water_inputdir,params.flooding_volume)
        varname = params.flooding_volume_varname  
    elif (temp_distrib == 'runoff'):
        pattern_file = os.path.join(params.water_inputdir,params.pnet)
        varname = params.pnet_varname      
    elif (temp_distrib == 'surface_runoff'):
        pattern_file = os.path.join(params.water_inputdir,params.surface_runoff)
        varname = params.surface_runoff_varname       
    elif (temp_distrib == 'land_runoff'):
        pattern_file = os.path.join(params.water_inputdir,params.land_runoff)
        varname = params.land_runoff_varname        
    #elif (temp_distrib == 'square_sro'):
    #    pattern_file = os.path.join(params.water_inputdir,params.surface_runoff)
    #    varname = params.surface_runoff_varname
    #    square = True     
    elif (temp_distrib == 'square_precip'):
        pattern_file = os.path.join(params.water_inputdir,params.precipitation)
        varname = params.precipitation_varname
        square = True
    else:
        print("Unknown temporal pattern provided in calculate_disaggregated_load: " + temp_distrib)
        raise MyError("Only discharge, flooding_volume, surface_runoff and square_precip are available!")
        
    try:
        pattern = MFDataset(pattern_file, 'r')
        #print('pattern file ' + pattern_file + ' is open!')
        # Convert dates of the pattern NETCDF data to decimal years
        ncyears = manip_netcdf.convert_numdate2year(pattern.variables['time'][:], pattern.variables['time'].units)
        #print('Dates in pattern file:')
        #print(str(ncyears))
        ibegin,iend = general_func.find_within_range(ncyears, year_start, year_end)

        # Apply pattern to yearly time series
        year_prev = year_start
        iyear = 0                    #index of the current year in the yearly time series
        load = load_yearly[iyear][1] #current yearly load
        #print('load  ' + str(load))
        year_dates = []              #time steps available for the current year
        year_vals = []               #pattern variable values for the current year
    
        for idate in range(ibegin, iend):
            year_dec = ncyears[idate]
            #print('year_dec  ' + str(year_dec) + '  ' + str(int(year_dec)))
            #print('year_start  ' + str(year_start))
            #print('year_prev  ' + str(year_prev))

            # If the date is inferior to the 1st simulated year, pass
            if (int(year_dec) < year_start):
                continue

            # Else start calculating pattern
            else:
                # If we are still in the same year, continue calculating yearly pattern
                if (int(year_dec) == year_prev):
                    year_dates.append(year_dec)
                    # Added to ignore negative values (ex. surface runoff): change?
                    pattern_val = max(0., pattern.variables[varname][idate, ipointer.get_lat_index(), ipointer.get_lon_index()])
                    if not(square):
                        #year_vals.append(pattern.variables[varname][idate, ipointer.get_lat_index(), ipointer.get_lon_index()])
                        year_vals.append(pattern_val)
                    else:
                        #year_vals.append(pattern.variables[varname][idate, ipointer.get_lat_index(), ipointer.get_lon_index()]**2)
                        year_vals.append(pattern_val**2)

                # Else apply pattern for previous year and start one for new year
                else:
                    for item in range(len(year_dates)):
                        if (sum(year_vals) > 0.0):
                            val = load * len(year_dates) * year_vals[item] / float(sum(year_vals))
                        else: #add warning if the pattern is equal to zero the whole year?
                            val = 0.0
                        load_dis.append([year_dates[item], val])

                    iyear += 1
                    load = load_yearly[iyear][1]
                    year_dates = [year_dec]
                    pattern_val = pattern.variables[varname][idate, ipointer.get_lat_index(), ipointer.get_lon_index()]
                    if not(square):
                        #year_vals = [pattern.variables[varname][idate, ipointer.get_lat_index(), ipointer.get_lon_index()]]
                        year_vals = [pattern_val]
                    else:
                        #year_vals = [pattern.variables[varname][idate, ipointer.get_lat_index(), ipointer.get_lon_index()]**2]
                        year_vals = [pattern_val**2]

                year_prev = int(year_dec)

        # Apply pattern to last year
        #print('year_dates  ' + str(year_dates))
        #print('year_vals  ' + str(year_vals))
        for item in range(len(year_dates)):
            if (sum(year_vals) > 0.0):
                val = load * len(year_dates) * year_vals[item] / float(sum(year_vals))
            else: #add warning if the pattern is equal to zero the whole year?
                val = 0.0
            load_dis.append([year_dates[item], val])

        pattern.close()
    
        return load_dis

    except:
        # If no temporal pattern file is provided, returns the yearly interpolated time series
        # We probably don't want that to happen... Apply 'uniform' pattern?
        print("No temporal pattern information is available in " + pattern_file)
        print("Yearly time series is returned.")
        return load_yearly

    
def get_dynamic_gridinfo(params,pointer1,filename,varname,temp_distrib=None,outputfile=None):
    '''
    Reads one variable (varname) for the whole simulation period and puts its values in a time dependent object.
    Reads from a NETCDF file, extracts the right time steps and applies a temporal pattern (to describe seasonnality) if asked for.
    Returned is dict object with cell number as key and a list of [time, value]. 
    '''

    obj_out = {}
    # Make empty list for all cells.
    for item in range(len(pointer1)):
        icell = pointer1[item].get_local_index()
        obj_out[icell] = []

    # Try to open dataset and apply pattern if asked for
    if True:
        if not "*" in filename:
          ncdata = Dataset(filename, 'r')
        elif '*' in filename:

          searchstr = copy.deepcopy(os.path.basename(filename)).replace('.nc', '')+'.zip'
          zip_list = directory.get_files_with_str(os.path.dirname(filename), searchstr)
          for zip_file in zip_list:
            zipfilePath = (zip_file)
            zip = zipfile.ZipFile(zipfilePath)
            zip.extractall(os.path.dirname(filename))
            zip.close()
          print(filename)
          ncdata = MFDataset(filename)
          searchstr = copy.deepcopy(os.path.basename(filename)).replace('.zip', '.nc')
          nc_list = directory.get_files_with_str(os.path.dirname(filename), searchstr, exclude=['*.zip*'])
         
          print(nc_list)
          for nc_file in nc_list:
            os.remove(nc_file)
        # Determine the years that are given for this parameter
        nctimes = ncdata.variables['time']
        years = manip_netcdf.convert_numdate2year(nctimes[:], nctimes.units)
        #print('years  ' + str(years))
    
        # Find begin and end point for this time period
        istart,iend = general_func.find_within_range(years,params.starttime,params.endtime)
        #print(varname)
        #print("istart,iend")
        #print(istart,iend)
        #print(min(iend+1,len(nctimes[:])))
        all_dat_period = ncdata.variables[varname][istart:min(iend+1,len(nctimes[:])),:,:] #wj
        
        # For each cell, create temp time series
        for item in range(len(pointer1)):
            icell = pointer1[item].get_local_index() # LV 12-07-2017
            out_cell = []
            datlist = all_dat_period[:, pointer1[item].get_lat_index(),pointer1[item].get_lon_index()] #wj
            #for idate in range(istart,iend):
            for idate in range(istart,min(iend+1,len(nctimes[:]))): #LV 20-02-2018
                #out_cell.append([years[idate],
                #                 ncdata.variables[varname][idate, pointer1[item].get_lat_index(), pointer1[item].get_lon_index()]])
                out_cell.append([years[idate], datlist[idate-istart]]) #wj
            #print('out_cell for item ' + str(item) + '  ' + str(out_cell))

            # If no pattern asked for, directly add time series to obj_out
            if (temp_distrib == None):
                obj_out[icell] = out_cell

            # Else create yearly time series and apply temporal pattern
            else:
                #out_cell_yearly = calculate_yearly_load(out_cell)
                out_cell_yearly = calculate_yearly_load(out_cell, int(params.starttime), int(params.endtime))
                out_cell_pattern = calculate_disaggregated_load(params, out_cell_yearly, pointer1[item], temp_distrib)
                obj_out[icell] = out_cell_pattern
        
    # If does not work, return the object with time_start en time_end with the value 0.0
    else:
        for item in range(len(pointer1)):
            icell = pointer1[item].get_local_index()
            obj_out[icell].append([params.starttime,0.0]) #add default value to function's arguments?
            obj_out[icell].append([params.endtime,0.0])
        print("No grid information is used for filename " + filename)
        print("Default value for the whole simulation period: 0.0")          
        
    if (outputfile == None):
        # Return the values
        #print('\nReturned object:')
        #print(str(obj_out))
        return obj_out
    else:
        # Dump obj_out to outputfile with pickle
        fp = open(outputfile,"wb")
        pickle.dump(obj_out,fp,-1)
        fp.close()


############
### Test ###
############
if __name__ == "__main__":

    import general_path

    try:
        # Test find within range
        years_test = [1887,1899,1900,1901,1903,1905,1907,1910]
        print("years_test  " + str(years_test))
        print("Test find_within_range in years test from 1900 to 1905:")
        print((general_func.find_within_range(years_test, 1900, 1905)))
        
        # Create a list of pointer objects
        pointer1_tmp = []
        npointers = 8
        nlon = 4
        for item in range(npointers):
            ilat = int(item/nlon)
            ilon = item - (ilat*nlon)
            pointer1_tmp.append(pointer_class.Pointer(index_cell=item,local_index=item,ilat=ilat,ilon=ilon,iorder=item))
            
        pointer1 = sorted(pointer1_tmp, key=lambda pointer_obj: pointer_obj.iorder)

        # Create test parameters
        params = Input_test()
        params.starttime = 1901
        params.endtime = 1901
        params.water_inputdir = '.'
        #params.discharge = 'pattern_test.nc'
        params.discharge = 'truc.nc'
        params.discharge_varname = 'qcmon'

        # Get dynamic grid info
        out_pattern = get_dynamic_gridinfo(params,pointer1,'test.nc','vartest',temp_distrib='discharge',outputfile=None)
        #out_pattern = get_dynamic_gridinfo_netcdf(params,pointer1,'truc.nc','vartest',temp_distrib='discharge',outputfile=None)

    except MyError as val:
        val.write()
    except (optparse.OptionValueError,
            Exception) as val:
        print(str(val))
        print(sys.exc_info()[0])
        print(traceback.print_exc())
    except:
        print('***** ERROR ******')
        print('get_dynamic_gridinfo.py failed.')
        print(sys.exc_info()[0])
        print(traceback.print_exc())
