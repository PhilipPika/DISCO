# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/retention_all.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os
import ascraster
import aggregate
import accuflux
import math


def temperature_factor(temp):
    '''
    Relation between temperature and retention (in rivers).
    With a temperature of 20.0 degrees Celcius, 1.0 is returned.
    The base number 1.06 (Marce and Armengol, 2009)
    We take air temperature as an approximation for water temperature.
    So air temperature has only not negative value
    '''
    return (1.06)**(max(0.0,temp)-20.0)

def retention_lake(params,tau):

    # Use model 4 from Brett and Benjamin (2008)
    # Tpout = Tpin/(1.0 + k4*(tau**x4)). Then retention_lake = (Tpin - Tpout)/Tpin = 1 - Tpout/Tpin = 1 - (1/(1.0 + k4*(tau**x4))
    # retention_lake = 1.0 - (1.0/(1.0 + k4*(tau**x4)))
    # with tau the residence time of the water in the lake, and k4 and x4 
    # calibrated factors.
    try:
        return 1.0 - (1.0/(1.0 + params.k4*(tau**params.x4))) 
    except ZeroDivisionError:
        return 0.0

def calculate_retention(params,hydraulic_load,water_body,temperature):

    # Calculate retention for a cell.
    # Retention is 1.0 - exp(-vf/Hl)
    if (water_body > 0.0):
        if (params.luse_temperature):
            temp_fac = temperature_factor(temperature)
        else:
            temp_fac = 1.0

        vf_temp_PO4 = params.P_vf * temp_fac
        try:
            return 1.0 - math.exp(- vf_temp_PO4 / hydraulic_load)
        except ZeroDivisionError:
            # when hydraulic load is zero, retention of the water body is equal to one.
            # So average over the cell is than flood_fraction *1.0 + (1-flood_fraction) * 0.0 = flood_fraction
            return 1.0
    else:
        return 0.0


def P_accuflux(params,mask,lake_icell_dict,lddmap,Pload,retention):

    '''
    ACCUFLUX
    Based on Gerard van Drecht, RIVM, Bilthoven, February 2002
    April 2002: accu fraction flux and accu fraction state added
    
    PURPOSE:
    Accumulate flow of mass, given by a material map [grid],
    along a flowpath, given by a local drain direction(LDD) map [grid].
    The output map contains accumulated flux [grid].
    All files are ascii-grid formatted.
    
    Lddmap: network of local drain directions
    flow_direction can be coded as follows:
    UNH/GRDC                       PCraster(LDD)
    32 64 128   meaning: NW N NE       7 8 9
    16  -   1             W -  E       4 5 6
     8  4   2            SW S SE       1 2 3
    We use the PCraster LDD-format; negative and zero values are assumed
    '''
    # Make output rasters
    accuPload = ascraster.duplicategrid(Pload)

    # Make order of all the cells.
    ordermap = accuflux.determine_order(lddmap)

    # Make a pointer for the calculation order
    pointer1 = []
    for icell in range(ordermap.length):
        pointer1.append([int(ordermap.get_data(icell)),icell])
    # Sort the pointer array so that lowest order is first calculated.
    pointer1.sort()

    # Loop over all cells:    
    for item in range(len(pointer1)):
        icell = pointer1[item][1]
        # Get mass flux at the cell.
        Pload_cell = accuPload.get_data(icell)
        if (Pload_cell == None):
            # No data value
            continue
        if Pload_cell < 0.0:
            print("Negative mass flow of " + str(Pload_cell) + " at cell: " + str(icell))           
            # Go to next cell.
            continue
        

        ret_cell = retention.get_data(icell,0.0)
        try:
            loutflow = (lake_icell_dict[icell].outcell == 1)
        except KeyError:
            loutflow = -1
        if (loutflow == 1): setattr(lake_icell_dict[icell],"load_in",Pload_cell)        
        Pload_cell = (1.0 - ret_cell) * Pload_cell
        if (loutflow == 1): setattr(lake_icell_dict[icell],"load_out",Pload_cell)
        accuPload.set_data(icell,Pload_cell)
        # Get direction of the ldd 
        direction = int(lddmap.get_data(icell,-1))
        if (direction < 1 or direction > 9 or direction == 5 ):
            pass
        else:
            icell_loop = accuflux.goto_nextcell(icell,direction,lddmap.ncols,mask=lddmap.mask)
            if (icell_loop > -1):
                prev = accuPload.get_data(icell_loop,0.0)
                accuPload.set_data(icell_loop,prev + Pload_cell)

    return accuPload


def calculate(params,mask,mouth_dict,lake_icell_dict,basin,discharge,lddmap,Pload):
    '''
    This function calculates the retention in the grid cells.
    '''

    # Read in the PCGLOBWB files
    water_storage_file = params.water_storage
    water_storage = ascraster.Asciigrid(ascii_file=water_storage_file,\
                                    mask=mask,numtype=float)

    flooding_depth_file = params.flooding_depth
    flooding_depth = ascraster.Asciigrid(ascii_file=flooding_depth_file,\
                                    mask=mask,numtype=float)  

    flooding_fraction = ascraster.Asciigrid(ascii_file=params.flooding_fraction,\
                                    mask=mask,numtype=float)

    fraction_water_grid = ascraster.Asciigrid(ascii_file=params.fraction_water,\
                                  mask=mask,numtype=float)          

    cellarea_grid = ascraster.Asciigrid(ascii_file=params.cellarea,\
                                    mask=mask,numtype=float)

    lakeid_grid = ascraster.Asciigrid(ascii_file=params.lakeid,\
                                    mask=mask,numtype=int)

    outlakeid_grid = ascraster.Asciigrid(ascii_file=params.outlakeid,\
                                    mask=mask,numtype=int)

    water_area_grid = ascraster.Asciigrid(ascii_file=params.water_area,\
                                    mask=mask,numtype=int)

    # Water volume flooding area
    volume_flooding_grid = ascraster.duplicategrid(flooding_depth)
    volume_flooding_grid.multiply(flooding_fraction,default_nodata_value =  0.0)
    volume_flooding_grid.multiply(cellarea_grid,default_nodata_value =  0.0)
    volume_flooding_grid.write_ascii_file(os.path.join(params.outputdir,"volume_flooding_grid.asc"))

    # Water volume of river.
    volume_river_grid = ascraster.duplicategrid(water_storage)
    volume_river_grid.write_ascii_file(os.path.join(params.outputdir,"volume_river_grid.asc"))

    # Area water of river.
    area_river_grid = ascraster.duplicategrid(cellarea_grid)
    area_river_grid.multiply(fraction_water_grid,default_nodata_value =  0.0)
    area_river_grid.write_ascii_file(os.path.join(params.outputdir,"area_river_grid.asc"))

    # Correct river volume and area when there is a reservoir or a lake.
    for icell in range(lakeid_grid.length):
        id = lakeid_grid.get_data(icell,-1)
        if (id > 0):
           # It is a lake or reservoir
           # Check whether it is a outlet of a lake or reservoir
           id = outlakeid_grid.get_data(icell,-1)
           if (id > 0):
               # It is an outlet. So put lake properties into the gridcell
               volume_river_grid.set_data(icell,water_storage.get_data(icell,0.0))
               area_river_grid.set_data(icell,water_area_grid.get_data(icell,0.0))
           else:
               # It is not an outlet. So remove the calculated river properties.
               volume_river_grid.set_data(icell,0.0)
               area_river_grid.set_data(icell,0.0)

    volume_river_grid.write_ascii_file(os.path.join(params.outputdir,"volume_river_grid_with_lakes.asc"))
    area_river_grid.write_ascii_file(os.path.join(params.outputdir,"area_river_grid_with_lakes.asc"))

    # Depth water of river in metres
    depth_river_grid = ascraster.duplicategrid(volume_river_grid)
    depth_river_grid.divide(area_river_grid,default_nodata_value =  0.0)
    depth_river_grid.write_ascii_file(os.path.join(params.outputdir,"depth_river_grid.asc"))

    # Change water_storage from m3 into km3
    volume_river_grid.multiply(1.0e-9)
    volume_flooding_grid.multiply(1.0e-9)

    # Residence time is the volume divided by discharge. 
    # Water storage is in km3. Discharge is km3/yr. So residence time is in years.
    residence_time_river_grid = ascraster.duplicategrid(volume_river_grid)
    residence_time_flooding_grid = ascraster.duplicategrid(volume_flooding_grid)
   
    # Calculate residence times by dividing water_storage through discharge.
    residence_time_river_grid.divide(discharge,default_nodata_value =  0.0)
    residence_time_flooding_grid.divide(discharge,default_nodata_value =  0.0,maximum=params.max_residence_time_per_cell)

    # If residence times are very high (e.g in dry areas),
    # set the residence time to a maximum
    # Apply maximum on gridcells which are not the outlet of a reservoir or lake.
    # In the flooding map there are no lakes/reservoirs
    for icell in range(residence_time_river_grid.length):
        resi = residence_time_river_grid.get_data(icell,-1)
        if (resi > params.max_residence_time_per_cell):
            # Check whether this is a outlet of a lake or reservoir
            try:
                outcell = lake_icell_dict[icell].outcell
                if (outcell == 1):
                    # Outlet of a lake/reservoir
                    if (resi > params.max_residence_time_per_lake):
                        resi = params.max_residence_time_per_lake
                else:
                    resi = params.max_residence_time_per_cell
                residence_time_river_grid.set_data(icell,resi)
            except KeyError:
                # This is no reservoir or lake
                residence_time_river_grid.set_data(icell,params.max_residence_time_per_cell)
    residence_time_river_grid.write_ascii_file(os.path.join(params.outputdir,"residence_time_river_grid.asc"))
    residence_time_flooding_grid.write_ascii_file(os.path.join(params.outputdir,"residence_time_flooding_grid.asc"))

    # Put residence time of the water body (lakes and reservoirs) in the grid.
    for icell in lake_icell_dict:
        if (lake_icell_dict[icell].outcell == 1):
            lake_icell_dict[icell].residence_time = residence_time_river_grid.get_data(icell,0.0)

    # Calculate the water volume in the grid cell
    volume_grid = ascraster.duplicategrid(volume_flooding_grid) 
    volume_grid.add(volume_river_grid,default_nodata_value =  0.0)
    volume_grid.write_ascii_file(os.path.join(params.outputdir,"volume_grid.asc"))

    # Calculate the average residence time with as weighing factors the volumes of the flooding - river
    residence_time_grid = ascraster.duplicategrid(volume_river_grid)
    residence_time_grid.multiply(residence_time_river_grid,default_nodata_value =  0.0)
    temp_grid = ascraster.duplicategrid(volume_flooding_grid)
    temp_grid.multiply(residence_time_flooding_grid,default_nodata_value =  0.0)
    residence_time_grid.add(temp_grid,default_nodata_value =  0.0)
    residence_time_grid.divide(volume_grid,default_nodata_value = 0.0)
    del temp_grid        
            
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

    # Calculate the hydraulic load (Hl = depth/residence_time)
    hydraulic_load_flooding_grid = ascraster.duplicategrid(flooding_depth)
    hydraulic_load_flooding_grid.divide(residence_time_flooding_grid,default_nodata_value =  0.0)
    hydraulic_load_river_grid = ascraster.duplicategrid(depth_river_grid)
    hydraulic_load_river_grid.divide(residence_time_river_grid,default_nodata_value = 0.0)
    hydraulic_load_river_grid.write_ascii_file(os.path.join(params.outputdir,"hydraulic_load_river.asc"))
    hydraulic_load_flooding_grid.write_ascii_file(os.path.join(params.outputdir,"hydraulic_load_flooding.asc"))
    # Calculate the average hydraulic load with as weighing factors the volumes of the flooding - river
    hydraulic_load_grid = ascraster.duplicategrid(volume_river_grid)
    hydraulic_load_grid.multiply(hydraulic_load_river_grid,default_nodata_value =  0.0) 
    temp_grid = ascraster.duplicategrid(volume_flooding_grid)
    temp_grid.multiply(hydraulic_load_flooding_grid,default_nodata_value =  0.0)
    hydraulic_load_grid.add(temp_grid,default_nodata_value =  0.0)
    hydraulic_load_grid.divide(volume_grid,default_nodata_value =  0.0)
    del temp_grid 

    # Write hydraulic load to output file:
    hydraulic_load_grid.write_ascii_file(os.path.join(params.outputdir,"hydraulic_load.asc")) 

    temperature_grid = ascraster.Asciigrid(ascii_file=params.water_temperature,\
                                      mask=mask,numtype=float) 

    # Read endoreic lakes
    endo_lakes_grid = ascraster.Asciigrid(ascii_file=params.endo_lakes,\
                                      mask=mask,numtype=int) 

    # Calculate the retention for rivers.
    retention_river_grid = ascraster.duplicategrid(volume_river_grid)
    for icell in range(retention_river_grid.length):
        id = endo_lakes_grid.get_data(icell,0)
        if (id > 0):
            # This is an endoreic lake. So the calculation
            retention_river_grid.set_data(icell,1.0)
            continue

        vol = volume_river_grid.get_data(icell,-1)
        if (vol > 0):
            temperature = temperature_grid.get_data(icell,20.0)
            hydraulic_load = hydraulic_load_river_grid.get_data(icell,0.0)
            ret = calculate_retention(params,hydraulic_load,vol,temperature)
            retention_river_grid.set_data(icell,ret)
        else:
            retention_river_grid.set_data(icell,0.0)

    # Calculate the retention for flood planes.
    retention_flooding_grid = ascraster.duplicategrid(volume_flooding_grid)
    for icell in range(retention_flooding_grid.length):
        vol = volume_flooding_grid.get_data(icell,-1)
        if (vol > 0):
            temperature = temperature_grid.get_data(icell,20.0)
            hydraulic_load = hydraulic_load_flooding_grid.get_data(icell,0.0)
            ret = calculate_retention(params,hydraulic_load,vol,temperature)
            retention_flooding_grid.set_data(icell,ret)
        else:
            retention_flooding_grid.set_data(icell,0.0)

    # Calculate the average retention with as weighing factors the volumes of the flooding - river
    retention_grid = ascraster.duplicategrid(volume_river_grid)
    retention_grid.multiply(retention_river_grid,default_nodata_value =  0.0) 
    temp_grid = ascraster.duplicategrid(volume_flooding_grid)
    temp_grid.multiply(retention_flooding_grid,default_nodata_value =  0.0)
    retention_grid.add(temp_grid,default_nodata_value =  0.0)
    retention_grid.divide(water_storage,default_nodata_value =  0.0)
    del temp_grid 

    # Set the lake database
    for icell in range(outlakeid_grid.length):
        if (id > 0):
            # It is an outlet. So put lake properties into the gridcell
            # Recalculate the temperature
            temp = 0.0
            number_of_cells = len(lake_icell_dict[icell].cells)
            for item in range(number_of_cells):
                icell_lake = lake_icell_dict[icell].cells[item]
                temp      += temperature.get_data(icell_lake,0.0)
            # Take average over the cells.
            temperature_cell = temp/float(number_of_cells)
            ret_cell = retention_river_grid.get_data(icell,0.0)
            setattr(lake_icell_dict[icell],"retention",ret_cell)
            setattr(lake_icell_dict[icell],"fraction_water",-99)
            setattr(lake_icell_dict[icell],"temperature",temperature_cell)

    accuPload = P_accuflux(params,mask,lake_icell_dict,lddmap,Pload,retention_grid)

    return residence_time_grid,traveltime,retention_grid,accuPload,hydraulic_load_grid

