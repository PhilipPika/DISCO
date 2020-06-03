# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/residence_time.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os
import ascraster
import aggregate
import accuflux

def calculate(params,mask,mouth_dict,lake_icell_dict,basin,discharge,lddmap):
    '''
    This function calculates the residence time in the grid cells.
    '''

    # Read water storage of rivers and lakes/reservoirs.
    water_storage = ascraster.Asciigrid(ascii_file=params.water_storage,\
                                    mask=mask,numtype=float)

    # Read rivers and lakes/reservoirs areas as fraction of the cellarea.
    fraction_waterbody_grid = ascraster.Asciigrid(ascii_file=params.fraction_water,\
                                    mask=mask,numtype=float)          

    # Read water depth of flooding areas.
    flooding_depth = ascraster.Asciigrid(ascii_file=params.flooding_depth,\
                                    mask=mask,numtype=float)  

    # Read flooding areas as fraction of the cellarea.
    flooding_fraction = ascraster.Asciigrid(ascii_file=params.flooding_fraction,\
                                    mask=mask,numtype=float)

    # Read the cellarea
    cellarea_grid = ascraster.Asciigrid(ascii_file=params.cellarea,\
                                    mask=mask,numtype=float)

    # Read lake/reservoir id
    lakeid_grid = ascraster.Asciigrid(ascii_file=params.lakeid,\
                                    mask=mask,numtype=int)

    # Read outlet cell of lake/reservoir id. Only outlet cell is filled.
    outlakeid_grid = ascraster.Asciigrid(ascii_file=params.outlakeid,\
                                    mask=mask,numtype=int)

    # Read area of lakes/reservoirs.
    waterbody_area_grid = ascraster.Asciigrid(ascii_file=params.water_area,\
                                    mask=mask,numtype=int)

    # Water volume flooding area in m3
    volume_flooding_grid = ascraster.duplicategrid(flooding_depth)
    volume_flooding_grid.multiply(flooding_fraction,default_nodata_value =  -1.0)
    volume_flooding_grid.multiply(cellarea_grid,default_nodata_value =  -1.0)


    # Water volume of waterbody (river,lakes or reservoir) in m3 as left over of water_storage minus flooding volume.
    volume_waterbody_grid = ascraster.duplicategrid(water_storage)
    volume_waterbody_grid.substract(volume_flooding_grid,default_nodata_value =  -1.0)

    # Area water of river in m2.
    area_waterbody_grid = ascraster.duplicategrid(cellarea_grid)
    area_waterbody_grid.multiply(fraction_waterbody_grid,default_nodata_value =  -1.0)

    # Correct river volume and area when there is a reservoir or a lake.
    for icell in range(lakeid_grid.length):
        id = lakeid_grid.get_data(icell,-1)
        if (id > 0):
           # It is a lake or reservoir
           # Check whether it is a outlet of a lake or reservoir
           id_lake = outlakeid_grid.get_data(icell,-1)
           if (id_lake > 0):
               # It is an outlet. So put lake properties into the gridcell
               volume_waterbody_grid.set_data(icell,water_storage.get_data(icell,0.0))
               area_waterbody_grid.set_data(icell,waterbody_area_grid.get_data(icell,0.0))
           else:
               # It is not an outlet. So remove the calculated river properties.
               volume_waterbody_grid.set_data(icell,0.0)
               area_waterbody_grid.set_data(icell,0.0)
    
    #volume_waterbody_grid.write_ascii_file(os.path.join(params.outputdir,"volume_waterbody_grid_negative.asc"))

    # Correct all the volumes which are negative.
    for icell in range(volume_waterbody_grid.length):     
        if (volume_waterbody_grid.get_data(icell,1.0) < 0.0):
            #print "Negative waterbody volume: ", volume_waterbody_grid.get_data(icell),\
            #      flooding_depth.get_data(icell),flooding_fraction.get_data(icell),\
            #      cellarea_grid.get_data(icell),water_storage.get_data(icell),fraction_waterbody_grid.get_data(icell)
            # When this happenes, we assume that the flooding is not correct. So the waterbody volume is equal to the water storage.
            # Water volume flooding is zero.
            volume_waterbody_grid.set_data(icell,water_storage.get_data(icell,0.0))
            volume_flooding_grid.set_data(icell,0.0)

    # Write volume in m3 and area in m2 to outputfile
    volume_waterbody_grid.write_ascii_file(os.path.join(params.outputdir,"volume_waterbody_grid.asc"))
    volume_flooding_grid.write_ascii_file(os.path.join(params.outputdir,"volume_flooding_grid.asc"))
    area_waterbody_grid.write_ascii_file(os.path.join(params.outputdir,"area_waterbody_grid.asc"))

    # Depth water of river in metres
    depth_waterbody_grid = ascraster.duplicategrid(volume_waterbody_grid)
    depth_waterbody_grid.divide(area_waterbody_grid,default_nodata_value =  -1.0)
    depth_waterbody_grid.write_ascii_file(os.path.join(params.outputdir,"depth_waterbody_grid.asc"))

    # Calculate the water volume in m3 in the grid cell
    volume_grid = ascraster.duplicategrid(volume_flooding_grid)
    volume_grid.add(volume_waterbody_grid,default_nodata_value = -1.0)
    volume_grid.write_ascii_file(os.path.join(params.outputdir,"volume_grid.asc"))

    # Change water_storage from m3 into km3
    volume_waterbody_grid.multiply(1.0e-9)
    volume_flooding_grid.multiply(1.0e-9)
    volume_grid.multiply(1.0e-9)

    # Residence time is the volume divided by discharge. 
    # Calculate the residence time as volume water body divided by corrected discharge
    # Corrected discharge is the discharge minus the water volume of the flooded area.
    
    # If residence times are very high (e.g in dry areas),
    # set the residence time to a maximum
    # Apply maximum on gridcells which are not the outlet of a reservoir or lake.
    residence_time_grid = ascraster.duplicategrid(volume_flooding_grid)
    for icell in range(residence_time_grid.length):
        # Check whether it is a lake/reservoir
        id = lakeid_grid.get_data(icell,-1)
        if (id > 0):
            # It is a lake or reservoir
            # Check whether it is a outlet of a lake or reservoir
            id_lake = outlakeid_grid.get_data(icell,-1)
            if (id_lake > 0):
                try:
                    resi = lake_icell_dict[icell].residence_time
                    # It is an outlet. So put lake properties into the gridcell         
                    if (lake_icell_dict[icell].residence_time > params.max_residence_time_per_lake):
                        # Change also the dictionary value of this waterbody
                        lake_icell_dict[icell].residence_time = params.max_residence_time_per_lake
                        resi = params.max_residence_time_per_lake
                except KeyError:
                     print(id , " has no residence time")
                     resi = 0.0           
            else:
                # Lake or reservoir but no outflow cell
                resi = 0.0
            residence_time_grid.set_data(icell,resi)
        else:
            # This is no reservoir or lake
            dis = discharge.get_data(icell,0.0)
            volume_flooded = volume_flooding_grid.get_data(icell,0.0)
            volume_waterbody = volume_waterbody_grid.get_data(icell,0.0)
            discharge_corrected = dis - volume_flooded
            if (dis <= 0.0):
                resi = 0.0
            elif (discharge_corrected >  volume_waterbody or volume_flooded <= 0.0):
                # Here we assume that there is enough water to feed the channel.
                resi = volume_waterbody/discharge_corrected
            else:
                # There is not enough water to feed the channel.
                # Now we take an volume weighed average of both residence times.
                volume_total = volume_waterbody + volume_flooded
                if (volume_total > 0.0):
                    try:
                        frac_flooded = volume_flooded/volume_total
                        resi_flooded = frac_flooded * volume_flooded/dis
                    except ZeroDivisionError:
                        frac_flooded = 0.0
                        resi_flooded = 0.0
                    try:
                        frac_waterbody = volume_waterbody/volume_total
                        resi_waterbody = frac_waterbody * volume_waterbody/dis
                    except ZeroDivisionError:
                        frac_waterbody = 0.0
                        resi_waterbody = 0.0

                    resi = resi_flooded + resi_waterbody
                else:
                    resi = 0.0
            resi = min(resi,params.max_residence_time_per_cell)
            residence_time_grid.set_data(icell,resi)
        
          
    # Write residence_time to output file:
    residence_time_grid.write_ascii_file(os.path.join(params.outputdir,"residence_time.asc"))     
 
    # Calculate the travel time of each grid cell to the river mouth
    traveltime = accuflux.accutime(lddmap,residence_time_grid,lake_icell_dict)
    traveltime.write_ascii_file(os.path.join(params.outputdir,"traveltime.asc")) 

    # Aggregate traveltime over the river basin.
    aggregate.aggregate_grid(basin,traveltime,mouth_dict,item="sum_traveltime")

    # Calculate the average traveltime for each river basin
    for key in list(mouth_dict.keys()):
        try:
            mouth_dict[key].avg_traveltime = mouth_dict[key].sum_traveltime/mouth_dict[key].number_of_cells 
        except ZeroDivisionError:
            mouth_dict[key].avg_traveltime = 0.0
        except AttributeError:
            pass

    return residence_time_grid,traveltime
