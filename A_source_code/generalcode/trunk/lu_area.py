# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/lu_area.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os
import ascraster
import aggregate

    
def calculate(params,mask,mouth_dict,basin):
    '''
    Calculates the area of each grass and arable land use cell.
    It returns two ascraster object with the area grass and area arable land in ha
    '''
       
    # Read grass files
    gras1 = ascraster.Asciigrid(ascii_file=params.pgrass,mask=mask,numtype=float)
    if params.ldebug: print(params.pgrass  + " has been read.")

    gras2 = ascraster.Asciigrid(ascii_file=params.igrass,mask=mask,numtype=float)
    if params.ldebug: print(params.igrass  + " has been read.")

    # Read crop files
    crop1 = ascraster.Asciigrid(ascii_file=params.cropmixp,mask=mask,numtype=float)
    if params.ldebug: print(params.cropmixp  + " has been read.")

    crop2 = ascraster.Asciigrid(ascii_file=params.croppasp,mask=mask,numtype=float)
    if params.ldebug: print(params.croppasp  + " has been read.")

    # Read land area
    landarea = ascraster.Asciigrid(ascii_file=params.landarea,mask=mask,numtype=float)
    if params.ldebug: print(params.landarea  + " has been read.") 
    
    # Aggregate landarea on river basin scale
    aggregate.aggregate_grid(basin,landarea,mouth_dict,item="landarea")
    
    # Make a copy (alias not a deepcopy) of the return values
    area_grass = gras1
    area_crop  = crop1
    area_natural = landarea
    
    # Calculate the total grass and crop land in ha. Because the LU file are in percentage, and the 
    # landarea is in km2, the product of (crop1_cell+crop2_cell) * landarea_cell is then in ha.
    for icell in range(crop1.length):
        crop1_cell = crop1.get_data(icell,0.0)
        crop2_cell = crop2.get_data(icell,0.0)
        gras1_cell = gras1.get_data(icell,0.0)
        gras2_cell = gras2.get_data(icell,0.0)
        landarea_cell = landarea.get_data(icell,0.0)
        # Put result in crop1 and gras1
        crop_area = landarea_cell * (crop1_cell + crop2_cell)
        gras_area = landarea_cell * (gras1_cell + gras2_cell)
        if ((crop_area + gras_area) - (100.0 * landarea_cell) > 1.):
            print("Crop and grass area is too large: ",crop_area + gras_area, 100. * landarea_cell, (crop_area + gras_area) - (100.0 * landarea_cell))
            try:
                fac = (100* landarea_cell)/(crop_area + gras_area)
            except ZeroDivisionError:
                fac = 0.0
            crop_area *= fac
            gras_area *= fac
            
        nat_area = max(0.0,(100.0 * landarea_cell) - crop_area - gras_area)
        # Put result in the maps in km2 (factor 0.01)
        area_crop.set_data(icell,0.01 * crop_area)
        area_grass.set_data(icell,0.01 * gras_area)
        area_natural.set_data(icell,0.01 * nat_area)
    
    # Aggregate landarea on river basin scale
    aggregate.aggregate_grid(basin,area_crop,mouth_dict,item="croparea")
    aggregate.aggregate_grid(basin,area_grass,mouth_dict,item="grasarea")
    aggregate.aggregate_grid(basin,area_natural,mouth_dict,item="natarea")

    # Calculate the number of cells in a region/basin
    one = ascraster.duplicategrid(landarea)
    # Set all cells to one.
    one.add_values(one.length*[1.0])
    aggregate.aggregate_grid(basin,one,mouth_dict,item="number_of_cells")

    # Delete the grid objects which are used anymore
    del crop1
    del crop2
    del gras1
    del gras2
    del landarea
    del one
   
    return area_grass,area_crop,area_natural
            
