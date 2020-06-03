# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/check_nodata.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************


# Import own general modules
from error import *
import ascraster

def check_nodata(filename,mask):
    '''
    This functions tests whether there are no nodata values given in this grid.
    When nodata is found, it will raise an error.
    '''

    # Open grid file
    grid = ascraster.Asciigrid(ascii_file=filename,mask=mask)

    # Check whether there is nodata.
    for icell in range(grid.length):
        val = grid.get_data(icell)
        if (val == None):
            # This should not happpen!
            # This should be corrected.
            raise MyError(filename + " has nodata values which are not allowed.")

