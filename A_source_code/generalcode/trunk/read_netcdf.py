# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/read_netcdf.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Reads information from netcdf file
# and writes the information to ascii grid output file

import sys
import os
import traceback

# Import Global ("generalcode") modules
# Generic code files ideally should be placed in the directory "generalcode"
# subfolder, but they may also be placed at the same base folder
__general = '/home/arthur/globalnutrients/generalcode/trunk'
#__general = os.path.join(r"../../..", 'generalcode/trunk')
if os.path.exists(__general): 
    sys.path.insert(0, __general)
    print(__general + " is added to the python search path for modules.") 

# General python modules
try:
    import numpy
    from netCDF4 import Dataset as DS
    from netcdftime import utime
    from datetime import  datetime
except ImportError:
    print("One of the general python modules is not installed.")
    sys.exit(2)

# Import generalcode modules.
import ascraster
import my_sys

def read_netcdf(netcdf_filename,varname,outputdir,outfile,year):
    '''
    Read out of data "varname" out of netcdf file.
    If year==None, then all years are written to outputdirectory. 
    '''

    # Location and name specification of netcdf files.
    #pcglobwb_dir = "PCR_Routing/netCDF"
    #basename_netcdf = "pcrglobwb_CRU21_30min_"
    #extension_netcdf = "_yearly.nc"
    #outputdir = "pcr_routing_out"

    # Get the command line arguments
    #year = sys.argv[1]
    #varname = sys.argv[2]
    #outfile = sys.argv[3]
    #extensie =sys.argv[4]

    #extension_netcdf = extensie + extension_netcdf

    if (not os.path.isdir(outputdir)):
        os.mkdir(outputdir)

    #netcdf_file = os.path.join(pcglobwb_dir,basename_netcdf + varname + extension_netcdf)
    output_file = os.path.join(os.path.join(outputdir,str(year)),outfile)

    if (not os.path.isfile(netcdf_filename)):
        print("***** ERROR *****")
        print("File is not found.")
        print(("Looked for file: " + netcdf_filename))
        print("***** ERROR *****")   
        raise IOError("File " + netcdf_filename + " is not found")

    # Load netcdf file
    nc=DS(netcdf_filename)

    # Check whether varname is in the netcdf file
    if (not varname in list(nc.variables.keys())):
        # There is a problem, varname is not available.
        print((varname + " is not found in file: " + netcdf_filename))
        print(("Parameters available: " + ",".join(map(str,list(nc.variables.keys())))))
        raise IOError("Parameter " + varname + " is not found in file " + netcdf_filename)

    # Get the time settings of the nc file
    try:
        cdftime = utime(nc.variables['time'].units,nc.variables['time'].calendar)
        start_year = cdftime.num2date(nc.variables['time'][0]).year
        end_year = cdftime.num2date(nc.variables['time'][-1]).year
    except AttributeError:
        # File is not time dependent.
        # Close netcdf file
        nc.close()
        read_netcdf_constant(netcdf_filename,varname,outputdir,outfile)


    if (year < start_year or year > end_year):
        print("***** ERROR *****")
        print(("Year: " + str(year) + " is not found in netcdf file " + netcdf_filename))
        print(("Years can be given between " + str(start_year) + " and " + str(end_year))) 
        print("***** ERROR *****")
        raise IOError("Year " + str(year) + " is outside the year range in file " + netcdf_filename)
    else:
        # Find the year in the netcdf file
        for i in range(len(nc.variables['time'])):
            if (cdftime.num2date(nc.variables['time'][i]).year == int(year)):
                iyear = i
                break
        else:
            print("***** ERROR *****")
            print(("Specified year " + str(year) + " is not found in netcdf file " + netcdf_filename))
            print("Years available: ")
            for i in range(len(cdftime.num2date(nc.variables['time']))):
                print(cdftime.num2date(nc.variables['time'][i]).year)
            raise IOError("Year " + str(int(year)) + " is not found in file " + netcdf_filename)

    # Get the variable.
    nc_var = nc.variables[varname][iyear,:,:]
    lat = nc.variables["latitude"][:]
    lon = nc.variables["longitude"][:]
    nrows = len(lat)
    ncols = len(lon)

    # Close netcdf file
    nc.close()

    dim=nc_var.shape

    if (dim[0] != nrows or dim[1] != ncols):
        print(("File: " + netcdf_filename + " is not consistent."))
        print(("In parameter file: nrows = ",nrows, " and ncols = ",ncols))
        print(("In netcdf file: nrows = ",dim[0], " and ncols = ",dim[1]))

    # General settings for ascii grid format  
    ncols = 720
    nrows = 360
    xll = -180
    yll = -90
    cellsize = 0.5
    nodata_value = -9999

    
    # Create a raster without a mask and filled with nodata.
    out_grid = ascraster.Asciigrid(ncols=ncols,nrows=nrows,\
                               xllcorner=xll, \
                               yllcorner = yll, \
                               cellsize = cellsize,\
                               nodata_value = nodata_value,\
                               numtype=float)            

    # Put the nc_var data into an ascraster object
    for i in range(nrows):
        for j in range(ncols):
            # Check whether there is no data
            if (type(nc_var[i,j]) != numpy.ma.core.MaskedConstant):
                icell = out_grid.get_index_from_coordin(lon[j],lat[i])
                out_grid.set_data(icell, nc_var[i,j])

    if (not os.path.isdir(os.path.join(outputdir,year))):
        os.mkdir(os.path.join(outputdir,year))
            
    out_grid.write_ascii_file(output_file)          

def read_netcdf_constant(netcdf_filename,varname,outputdir,outfile):
    '''
    Read out of data "varname" out of netcdf file.
    If year==None, then all years are written to outputdirectory. 
    '''

    if (not os.path.isdir(outputdir)):
        os.mkdir(outputdir)

    # Make output file
    output_file = os.path.join(outputdir,outfile)

    if (not os.path.isfile(netcdf_filename)):
        print("***** ERROR *****")
        print("File is not found.")
        print(("Looked for file: " + netcdf_filename))
        print("***** ERROR *****")   
        raise IOError("File " + netcdf_filename + " is not found")

    # Load netcdf file
    nc=DS(netcdf_filename)

    # Check whether varname is in the netcdf file
    if (not varname in list(nc.variables.keys())):
        # There is a problem, varname is not available.
        print((varname + " is not found in file: " + netcdf_filename))
        print(("Parameters available: " + ",".join(map(str,list(nc.variables.keys())))))
        raise IOError("Parameter " + varname + " is not found in file " + netcdf_filename)


    # Get the variable.
    nc_var = nc.variables[varname][:,:]
    try:
        lat = nc.variables["latitude"][:]
    except KeyError:
        lat = nc.variables["lat"][:]
    try:
        lon = nc.variables["longitude"][:]
    except KeyError:
        lon = nc.variables["lon"][:]
    nrows = len(lat)
    ncols = len(lon)

    # Close netcdf file
    nc.close()

    dim=nc_var.shape

    if (dim[0] != nrows or dim[1] != ncols):
        print(("File: " + netcdf_filename + " is not consistent."))
        print(("In parameter file: nrows = ",nrows, " and ncols = ",ncols))
        print(("In netcdf file: nrows = ",dim[0], " and ncols = ",dim[1]))

    # General settings for ascii grid format  
    ncols = 720
    nrows = 360
    xll = -180
    yll = -90
    cellsize = 0.5
    nodata_value = -9999

    # Create a raster without a mask and filled with nodata.
    out_grid = ascraster.Asciigrid(ncols=ncols,nrows=nrows,\
                               xllcorner=xll, \
                               yllcorner = yll, \
                               cellsize = cellsize,\
                               nodata_value = nodata_value,\
                               numtype=float)            

    # Put the nc_var data into an ascraster object
    for i in range(nrows):
        for j in range(ncols):
            # Check whether there is no data
            if (type(nc_var[i,j]) != numpy.ma.core.MaskedConstant):
                icell = out_grid.get_index_from_coordin(lon[j],lat[i])
                out_grid.set_data(icell, nc_var[i,j])

    if (not os.path.isdir(os.path.join(outputdir,year))):
        os.mkdir(os.path.join(outputdir,year))
            
    out_grid.write_ascii_file(output_file)          



if __name__ == '__main__':
   netcdf_filename = "/home/arthur/globalnutrients/fix_input/Pdeposition/depannual.nc"
   outputdir = "/home/arthur/out_netcdf"
   varname = "currentdep"
   outfile = varname
   year = 2000
               
   read_netcdf(netcdf_filename,varname,outputdir,outfile,year)

