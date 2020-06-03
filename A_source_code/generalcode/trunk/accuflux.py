# ******************************************************
## Revision "$LastChangedDate: 2018-11-15 16:34:39 +0100 (do, 15 nov 2018) $"
## Date "$LastChangedRevision: 11 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/accuflux.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import ascraster
import os
import sys

icold = [0,-1,0,1,-1,0,1,-1,0,1]
irowd = [0,1,1,1,0,0,0,-1,-1,-1]

def get_row_col_from_index(ind,ncols):
    '''
    Returns irow and icol of an index
    '''
    if (type(ncols) != int or type(ind) != int):
        raise MyError("Function get_row_col_from_index must have integers as argument")
    # Note, here a division betwen integers is needed.    
    irow = ind//ncols
    icol = ind - irow * ncols
    return irow,icol

def get_index_from_row_col(irow, icol, ncols):
    ind = irow * int(ncols) + icol
    return ind

def get_neighbours(id, ncols, mask=None):
    nbs = list()
    if mask==None:
        irow,icol = get_row_col_from_index(id,ncols)
    else:
        irow, icol = get_row_col_from_index(mask[id],ncols)
    for row in range(irow-1, irow+2):
        for col in range(icol-1, icol+2):
            if not (col==icol and row==irow):
                nbs.append(get_index_from_row_col(row, col, ncols))
    return nbs

def check_cell(id_new,id,ncols,mask=None):
    '''
    Check whether cell id1 and id2 are neighbours.
    It is possible that it returns a cell id, that falls outside the 
    raster (larger than length of raster). We need ncols for this.
    '''
    # Get column number and row number of these cells.
    if (mask == None):
        irow,icol = get_row_col_from_index(id,ncols)
        irow_new,icol_new = get_row_col_from_index(id_new,ncols)
    else:
        irow,icol = get_row_col_from_index(mask[id],ncols)
        irow_new,icol_new = get_row_col_from_index(mask[id_new],ncols)
    
    #print "COLS:",irow,icol,irow_new,icol_new
    if (icol - icol_new > 1 or icol_new - icol > 1):
        # This is not a neighbour
        # print "NO: col"
        return -1
    if (irow -irow_new > 1 or irow_new - irow > 1):
        # This is not a neighbour
        # print "NO: row"
        return -1
    if (icol_new < 0 or icol_new >= ncols):
        # This is not a neighbour
        # print "NO: outside range column"
        return -1
    if (irow_new < 0):
        # This is not a neighbour
        # print "NO: negative rownumber"
        return -1            
    # print "YES:"
    return id_new
        
def goto_nextcell(icell,direction,ncols,mask=None):
    '''
    Return the cell number of the first cell which is found in the direction of direction.
    When the grid is mask, this is difficult.
    When the grid is not masked, it is easy.
    Direction numbers:
        UNH/GRDC                       PCraster(LDD)
        32 64 128   meaning: NW N NE       7 8 9
        16  -   1             W -  E       4 5 6
         8  4   2            SW S SE       1 2 3
    We use the PCraster LDD-format; 
    '''

    if (direction == 5):
        # This is a pit. The neighbour of this cell is itself.
        return icell

    if (mask == None):
        icell_new = icell + irowd[direction]*ncols + icold[direction]
        id = check_cell(icell_new,icell,ncols)
        return id 
    else:
        # The grid is masked
        # icell_new is the number of the gridcell in the unmasked grid.
        icell_new = mask[icell] + irowd[direction]*ncols + icold[direction]
        try:
            if (icell_new > mask[icell]):
                # Find value icell_new in mask[icell:]
                id = mask.index(icell_new,icell+1,icell+ncols+2)
            else:
                icell_new_pointer = icell + irowd[direction]*ncols + icold[direction]
                if (icell_new_pointer > mask[-1]): icell_new_pointer = mask[-1]
                if (icell_new_pointer < 0): icell_new_pointer = 0

                # Find value icell_new in mask[icell_new_pointer:]
                id = mask.index(icell_new,icell_new_pointer,icell)
                
            # Check whether the id is a neighbour.
            id = check_cell(id,icell,ncols,mask=mask)
            
        except ValueError:        
            #print "The next cell is not in the mask for cell: " + str(icell) + " and direction " +str(direction)
            id = -1

        return id


def factor(icell,retentionmap):
    '''
    When retentionmap is given is will give the 1.0 - value of retentionmap in the cell icell.
    Otherwise it will return 1.
    '''
    if (retentionmap == None):
        return 1.0
    else:
        fac = retentionmap.get_data(icell,0.0)
        if (fac < 0.0):
            print("retentionmap is not correct. It contains negative value of " + str(fac) + " in cell: "+str(icell))
            # Change retention map.
            retentionmap.set_data(icell,0.0)
            return 0.0
        elif (fac > 1.0):    
            print("retentionmap is not correct. It contains positive value (> 1.0) of " + str(fac) + " in cell: "+str(icell))
            # Change retention map.
            retentionmap.set_data(icell,1.0)
            return 1.0
        else:
            return (1.0 - fac)
            
def accuflux(lddmap,fluxmap,retentionmap=None,negative_flux_possible=0):
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

    if (negative_flux_possible == 1):
        # Call the other version of accuflux
        return accuflux_negative(lddmap,fluxmap,retentionmap=retentionmap,negative_flux_possible=1)

    # Make output rasters
    accufluxmap = ascraster.duplicategrid(fluxmap)
    # Make nodata on all cells.
    accufluxmap.add_values(accufluxmap.length*[accufluxmap.nodata_value])
           
    # Loop over all cells:    
    for icell in range(lddmap.length):
        # Get mass flux at the cell.
        flux = fluxmap.get_data(icell)
        if (flux == None):
            # No data value
            continue
        if (flux < 0.0):
            if (negative_flux_possible == 0):
                print("Negative mass flow of " + str(flux) + " at cell: " + str(icell))
                # Go to next cell.
                continue
            
        end_of_path=False
        icell_loop = icell
        while not end_of_path:
            # Get direction of the ldd 
            direction = int(lddmap.get_data(icell_loop,-1))
            if (direction < 1):
                #print "Something is wrong with the lddmap (< 1). No river mouth found for icell: ",icell, " with load: ", flux, " and direction: ",direction
                # Go to next cell
                end_of_path=True
            elif (direction == 5):
                # Mouth of river basin is found.
                prev = accufluxmap.get_data(icell_loop,0.0)
                accufluxmap.set_data(icell_loop,prev + flux * factor(icell_loop,retentionmap))
                # Go to next cell
                end_of_path=True
            elif (direction > 9):
                #print "Something is wrong with the lddmap (> 9). No river mouth found for icell: ",icell, " with load: ", flux, " and direction: ",direction
                # Go to next cell
                end_of_path=True
            else:
                prev = accufluxmap.get_data(icell_loop,0.0)
                flux = flux * factor(icell_loop,retentionmap)
                accufluxmap.set_data(icell_loop,prev + flux)
                
                # Go to next cell
                icell_loop = goto_nextcell(icell_loop,direction,lddmap.ncols,mask=lddmap.mask)
                if (icell_loop < 0):
                    # Already printed a message. Go to next grid cell.
                    end_of_path=True

    return accufluxmap

def accuflux_negative(lddmap,fluxmap,retentionmap=None,negative_flux_possible=1):
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
    The flux is allowed to be negative, but the accufluxmap is only allowed to be positive!
    '''
    
    # Make output rasters
    accufluxmap = ascraster.duplicategrid(fluxmap)
    # Make nodata on all cells.
    accufluxmap.nodata_value = -1
    accufluxmap.add_values(accufluxmap.length*[accufluxmap.nodata_value])

    # Make order of all the cells.
    ordermap = determine_order(lddmap)

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
        flux = fluxmap.get_data(icell)
        if (flux == None):
            # No data value
            continue

        # Set current accufluxmap
        prev = accufluxmap.get_data(icell,0.0)
        flux = max(0.0,(flux+prev) * factor(icell,retentionmap))
        accufluxmap.set_data(icell,flux)
        
       
        # Bring mass to the next cell
        direction = int(lddmap.get_data(icell,-1))
        if (direction == 5):
            continue
        icell_next = goto_nextcell(icell,direction,lddmap.ncols,mask=lddmap.mask)
        if (icell_next >  -1):
            prev_next = accufluxmap.get_data(icell_next,0.0)
            accufluxmap.set_data(icell_next,prev_next + flux)

    return accufluxmap



def determine_order(lddmap):
    # Make output raster with order of the cell
    ordermap = ascraster.duplicategrid(lddmap)
    # Fill all cells with zero for the order and one for the valuemap.
    ordermap.add_values(ordermap.length*[0])

    # Loop over all cells:    
    for icell in range(lddmap.length):         
        end_of_path=False
        icell_loop = icell
        icell_loop_prev = icell
        while not end_of_path:
            lfirst = (icell_loop_prev == icell_loop)
            # Get direction of the ldd 
            direction = int(lddmap.get_data(icell_loop,-1))
            neighbour = int(ordermap.get_data(icell_loop_prev,0))
            if (direction < 1):
                # Go to next cell
                end_of_path=True
            elif (direction == 5):
                # Mouth of river basin is found.
                prev = int(ordermap.get_data(icell_loop,0))
                if (lfirst):
                    if (prev < neighbour):
                        ordermap.set_data(icell_loop,neighbour + 1)
                else:
                    if (prev < neighbour + 1):
                        ordermap.set_data(icell_loop,neighbour + 1)
                # Go to next cell
                end_of_path=True
            elif (direction > 9):
                # Go to next cell
                end_of_path=True
            else:
                prev = int(ordermap.get_data(icell_loop,0))
                if (lfirst):
                    if (prev < neighbour):
                        ordermap.set_data(icell_loop,neighbour + 1)
                else:
                    if (prev < neighbour + 1):
                        ordermap.set_data(icell_loop,neighbour + 1)
                
                # Go to next cell
                icell_loop_prev = icell_loop
                icell_loop = goto_nextcell(icell_loop,direction,lddmap.ncols,mask=lddmap.mask)

                if (icell_loop < 0):
                    # Already printed a message. Go to next grid cell.
                    end_of_path=True

    return ordermap

def make_strahlerordermap(lddmap):
    '''
    '''
    cell = []
    source_array = [cell]*lddmap.length

    strahler = ascraster.duplicategrid(lddmap)
    order_map = determine_order(lddmap)
    max_order = order_map.maximum()
    order_map.write_ascii_file('order_map.asc')

    # Loop over all cells to create source list:	
    for o in range(0, max_order+1):
        ordermask = ascraster.create_mask('order_map.asc', o, 'EQ',numtype=int)

        for icell in ordermask: 	
            if (lddmap.get_data(icell) != None):
                if (o==0):
                    strahler.set_data(icell, 6)
                icell_next = goto_nextcell(icell,lddmap.get_data(icell),lddmap.ncols)
                if (source_array[icell_next]==[]):
                    source_array[icell_next] = [icell]
                else:
                    source_array[icell_next].append(icell)
                if (o!=0):
                    dum = list()
                    for source_icell in source_array[icell]:
                        dum.append(strahler.get_data(source_icell))
                    if (len(dum)>1):
                        if (all(x == dum[0] for x in dum)):
                            strahler.set_data(icell, dum[0]+1)
                        else:
                            strahler.set_data(icell, max(dum))
                    else:
                        strahler.set_data(icell, dum[0])
    strahler.write_ascii_file('real_strahler.asc')

def accutime(lddmap,timemap,lake_icell_dict):
    '''
    ACCUTIME
   
    PURPOSE:
    Accumulate residence time for each grid cell to the mouth of the river.

    The output map contains accumulated residence time [grid].
    When a endoreic lake is encountered then, residence time is to the lake and not to the sea.
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
    accutimemap = ascraster.duplicategrid(timemap)
          
    # Loop over all cells:    
    for icell in range(lddmap.length):
        # Get residence time at the cell.
        restime = timemap.get_data(icell)
        if (restime == None):
            # No data value
            continue
        if (restime < 0.0):
            print("Negative residence time of " + str(restime) + " at cell: " + str(icell))
            # Go to next cell.
            accutimemap.set_data(icell,0.0)
            continue
            
        end_of_path=False
        icell_loop = icell
        cumtime = restime
        while not end_of_path:
            # Get direction of the ldd 
            direction = int(lddmap.get_data(icell_loop,-1))
            # Check whether the cell is in a endoreic lake
            try:
                 lendo_lake = lake_icell_dict[icell].lendolake
            except KeyError:
                 lendo_lake = 0
            if (direction < 1 or direction > 9):
                #print "Something is wrong with the lddmap (< 1). No river mouth found for icell: ",icell, " with load: ", flux, " and direction: ",direction
                # Go to next cell
                end_of_path=True
            elif (direction == 5):
                # Mouth of river basin is found.
                end_of_path=True
            elif (lendo_lake == 1):
                # Endoreic lake found.
                end_of_path=True
            else:               
                # Go to next cell
                icell_loop = goto_nextcell(icell_loop,direction,lddmap.ncols,mask=lddmap.mask)
                if (icell_loop < 0):
                    # Already printed a message. Go to next grid cell.
                    end_of_path=True

            # Set accutime in grid or fetch new restime for next cell
            if (end_of_path):
                accutimemap.set_data(icell,cumtime)
            else:
                restime = timemap.get_data(icell_loop,0.0)
                cumtime += restime
       

    return accutimemap
