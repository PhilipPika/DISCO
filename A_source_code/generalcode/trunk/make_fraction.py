# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/make_fraction.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import ascraster

    
def make_fraction(params,source,total,divisionerror_outcome=0.0):
    '''
    Calculates the fraction of the source on the total for each grid cell.
    '''
       
    if params.ldebug: print("Start with calculation of make_fraction.")

    fraction = ascraster.duplicategrid(total)
    if (fraction.nodata_value == None):
        fraction.nodata_value = -1

    for icell in range(total.length):
        total_cell = total.get_data(icell,0.0)
        source_cell = source.get_data(icell,0.0)
        try:
            frac = source_cell/total_cell
        except ZeroDivisionError:   
            frac = divisionerror_outcome
            
        # Put result in raster
        fraction.set_data(icell,frac)

    return fraction
    
