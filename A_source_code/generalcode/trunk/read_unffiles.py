# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/read_unffiles.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
'''
read_unffile.py

A Python library containing functions to read IMAGE UNF files and store as ASCII raster files.

Initialy created on 13 april 2015
@author: beusena
'''
import os

import ascraster
from error import *
import numpy as np

def read_unffile(filename):
    '''
    Reads UNF files from IMAGE model and stores it as ascraster object.
    The example raster is used to convert the unformatted file to a ascraster.
    '''
    # Determine the datatype
    dt,number_of_grids = unf_datatype(filename)

    # Read unformatted file
    datavalues = np.fromfile(filename, dtype=dt).tolist()

    return datavalues

def read_unf_as_raster(filename,example_raster):
    '''
    Reads UNF files from IMAGE model and stores it as ascraster object.
    The example raster is used to convert the unformatted file to a ascraster.
    '''
    # Determine the datatype
    dt,number_of_grids = unf_datatype(filename)

    # Read unformatted file
    datavalues = np.fromfile(filename, dtype=dt).tolist()

    gridout = []
    for item in range(number_of_grids):
        gridout.append(ascraster.duplicategrid(example_raster))

    if (len(datavalues) == number_of_grids * gridout[0].length):
        # Size is okay
        if (number_of_grids == 1):
            gridout[0].values = copy.deepcopy(datavalues)
        else:
            # Datavalues are stored for each cell.
            for igrid in range(number_of_grids):
                for icell in range(gridout[0].length):
                    gridout[igrid].set_data(icell,datavalues[icell * number_of_grids + igrid])
    else:
        raise MyError("Size of example raster is not equal to unformatted file.",\
                      "Size of example raster: " + str(gridout.length),\
                      "Size of unformatted datafile: " + str(len(datavalues)) + " in file: " + filename)

    if (number_of_grids == 1):
        # Return the ascraster object (not as a list).
        return gridout[0]
    else:
        # Return all the grids as list
        return gridout

def unf_datatype(filename):
    '''
    Data type is determined on bases of the extension(s).
    Returned the type format for numpy and the number of classes in the unf file
    '''
    
    # Split filename into basename and extension.
    basename,extension = os.path.splitext(filename)

    # Find the data type in the extension of the filename
    if (extension.upper() == ".UNF0"):
        ctype = "f4"
    elif (extension.upper() == ".UNF1"):
        ctype = "i1"
    elif (extension.upper() == ".UNF2"):
        ctype = "i2"
    elif (extension.upper() == ".UNF4"):
        ctype = "i4"
    else:
        raise MyError("Unknown extension found for file: ",filename)

    # Check whether there are more than one value per cell stored.
    # When there are more than value stored, it will be behind the last . in the filename.
    fields = basename.split(".")
    try:
        number_of_grids = int(fields[-1])
        numtxt = fields[-1]
    except:
        numtxt = ""
        number_of_grids = 1
 
    # Combine the number of fields and the type
    return numtxt + ctype,number_of_grids


# Testing this module           
if (__name__ == "__main__"):
    
    print(unf_datatype("GFRAC_1970.19.UNF0"))
    print(unf_datatype("GAREALAND.UNF0"))
    print(unf_datatype("GR.UNF1"))
    print(unf_datatype("GCOUNTRY.UNF2"))
    print(unf_datatype("GLCT.UNF4"))
    print(unf_datatype("GLCT.unf2"))
    try:
        print(unf_datatype("GLCT.unf3"))
    except MyError as val:
        val.write()

    filename = r'/home/arthur/mestverdeling/in1900/GAREALAND.UNF0'
    dt,num = unf_datatype(filename)
    area = np.fromfile(filename, dtype=dt).tolist()
    print(area[:10])
    area = read_unffile(filename)
    print(area[:10])
 
