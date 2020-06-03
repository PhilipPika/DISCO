# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/hydraulic_load.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os
import ascraster

def calculate(params,mask,residence_time_grid,lake_icell_dict):
    '''
    This function calculates the hydraulic load in the water bodies.
    '''

    # Depth water of waterbody in metres (result of residence time function)
    depth_grid = ascraster.Asciigrid(ascii_file=os.path.join(params.outputdir,"depth_waterbody_grid.asc"),\
                                              mask=mask,numtype=float)   

    hydraulic_load_grid = ascraster.duplicategrid(residence_time_grid)
    hydraulic_load_grid.add_values(hydraulic_load_grid.length *[0.0])
    for icell in range(hydraulic_load_grid.length):
        depth = depth_grid.get_data(icell,0.0)
        resi_time = residence_time_grid.get_data(icell,0.0)
        try:
            Hl = depth/resi_time
        except ZeroDivisionError:
            Hl = 0.0
        hydraulic_load_grid.set_data(icell,Hl)

    # Write hydraulic load to output file:
    hydraulic_load_grid.write_ascii_file(os.path.join(params.outputdir,"hydraulic_load.asc"))

    return hydraulic_load_grid
