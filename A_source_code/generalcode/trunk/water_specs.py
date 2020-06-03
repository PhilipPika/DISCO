# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/water_specs.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os
import ascraster
import areacell

def calculate(params,mask,mouth_dict,basin,pnet,discharge):
    '''
    This function calculates the water area and the water body.
    Here the discharge of the mouth and qs is set and calculated.
    '''
    m2tokm2 = 1.0e-6
    mmtokm  = 1.0e-6
    kmtom   = 1.0e3

    # Read flooding areas as fraction of the cellarea.
    flooding_fraction = ascraster.Asciigrid(ascii_file=params.flooding_fraction,\
                                    mask=mask,numtype=float)

    # Read rivers and lakes/reservoirs areas as fraction of the cellarea.
    fraction_water_grid = ascraster.Asciigrid(ascii_file=params.fraction_water,\
                                  mask=mask,numtype=float)    

    # Read the cellarea
    cellarea_grid = ascraster.Asciigrid(ascii_file=params.cellarea,\
                                        mask=mask,numtype=float)
    cellarea_grid.multiply(m2tokm2)

    # Read lake/reservoir id
    lakeid_grid = ascraster.Asciigrid(ascii_file=params.lakeid,\
                                    mask=mask,numtype=int)

    # Calculate water area
    water_area = ascraster.duplicategrid(flooding_fraction)
    local_discharge = ascraster.duplicategrid(discharge)

    # Make nodata on all cells.
    water_area.add_values(water_area.length*[water_area.nodata_value])
    flooding_area = ascraster.duplicategrid(water_area)

    for icell in range(flooding_fraction.length):
        id = lakeid_grid.get_data(icell,-1)
        if (id > 0):
            waterbody_fraction_cell = fraction_water_grid.get_data(icell,0.0)
            flooding_fraction_cell  = 0.0
        else:
            flooding_fraction_cell = flooding_fraction.get_data(icell,0.0)
            waterbody_fraction_cell = fraction_water_grid.get_data(icell,0.0)
            waterbody_fraction_cell += flooding_fraction_cell
        if (waterbody_fraction_cell > 1.0):
            waterbody_fraction_cell = 1.0
        area = cellarea_grid.get_data(icell,0.0)
        water_area.set_data(icell,waterbody_fraction_cell * area)
        flooding_area.set_data(icell,flooding_fraction_cell * area)
        
        # Calculation of local discharge.
        discharge_cell = pnet.get_data(icell,0.0)
        if (discharge_cell > 0.):
            # Convert pnet [mm/yr] to km3 by multiplying with cellarea and 10e-6 (conversion from mm to km)
            flux = discharge_cell * area * mmtokm
            local_discharge.set_data(icell,flux)
        else:
            local_discharge.set_data(icell,0.0)

    # Take the discharge of the mouth of the river divided by the water area of the mouth cell. 
    for key in list(mouth_dict.keys()):
       mouth_dict[key].discharge_mouth = discharge.get_data(mouth_dict[key].index_grid,0.0)
       mouth_dict[key].waterarea = 0.0
       
    # Calculate the water body area for each river basin.
    # Create a start residence time with zero's.
    water_body = ascraster.duplicategrid(discharge)
    # Make nodata on all cells. As long as we can not calculate the waterarea, we use a waterbody of one instead of zero!
    water_body.add_values(water_body.length*[0])

    for icell in range(0,basin.length):
       basinid = basin.get_data(icell,-1)
       if (basinid > 0):
           try:
               area = water_area.get_data(icell,0.0)
               if (area > 0.0):                  
                   mouth_dict[int(basinid)].waterarea += area
                   water_body.set_data(icell,1)
                   mouth_dict[int(basinid)].flooding_area += flooding_area.get_data(icell,0.0)
           except KeyError:
               # We don't not do anything with this river.
               print("Basinid " + str(basinid) + " found in basinid map but not in river mouth file.")
    
    for key in list(mouth_dict.keys()):
       try:
           mouth_dict[key].qs = kmtom * mouth_dict[key].discharge_mouth/mouth_dict[key].waterarea
       except ZeroDivisionError:
            mouth_dict[key].qs  = 0.0

    # Write local discharge to output file
    local_discharge.write_ascii_file(os.path.join(params.outputdir,"discharge_local.asc"))

    return water_area, water_body
