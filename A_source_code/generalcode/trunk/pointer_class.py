# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/pointer_class.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import ascraster
import accuflux
from iround import *

class Pointer:
    '''
    Store the various pointer index of a gridcell in a class.
    '''
    def __init__(self,index_cell=None,local_index=None,ilat=None,ilon=None,iorder=None):
        # Index in the unmasked grid.
        self.index_cell = index_cell
        # Index in the masked grid.
        self.local_index = local_index
        # Row index or latitude index
        self.ilat = ilat
        # Column index or longitude index
        self.ilon = ilon
        # Order of the cell in the riverbasin. 
        # Higher number means closer to the mouth of the river.
        self.iorder = iorder

    def get_index(self):
        return self.index_cell

    def get_local_index(self):
        return self.local_index

    def get_lat_index(self):
        return self.ilat

    def get_lon_index(self):
        return self.ilon

    def get_iorder(self):
        return self.iorder


def make_pointer(params,rivmask=None):
    '''
    Make a list of Pointer items with all information from each cell.
    '''
    # Read ldd map 
    lddmap = ascraster.Asciigrid(ascii_file=params.ldd,mask = rivmask,numtype=int)

    # Calculate the order of each cell
    ordermap = accuflux.determine_order(lddmap)

    # Make a pointer for the calculation order
    pointer1 = []
    for icell in range(ordermap.length):
        pointer1.append([iround(ordermap.get_data(icell)),icell])
    # Sort the pointer array so that lowest order is first calculated.
    pointer1.sort()

    # Create a dummy grid which has no mask to get the lat-lon numbers.
    dummy = ascraster.Asciigrid(ncols = lddmap.ncols, nrows = lddmap.nrows, \
                                xllcorner = lddmap.xllcorner, yllcorner = lddmap.yllcorner,\
                                cellsize = lddmap.cellsize)

    # Put all in a list of Pointer objects
    pointer_class = []
    for item in range(len(pointer1)):
        local_index = pointer1[item][1]
        iorder =  pointer1[item][0]
        if (rivmask == None):
            index_cell = local_index
        else:
            index_cell = rivmask[local_index]
        ilat,ilon = dummy.get_row_col_from_index(index_cell)
        pointer_class.append(Pointer(index_cell=index_cell,local_index=local_index,\
                                     ilat=ilat,ilon=ilon,iorder=iorder))

    return pointer_class
