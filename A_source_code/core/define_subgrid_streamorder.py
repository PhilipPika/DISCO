# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/define_subgrid_streamorder.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# General python modules
import os

# Import own general modules
import ascraster
import general_class

def define_subgrid_streamorder_constant(params):
    '''
    This function defines the constant factors for the subgrid_streamorder (cell independent).
    The stream orders 1 to params.norder are defined, 
    where the last order is the first stream that is modelled by PCRGLOBWB.
    All output parameters are set in the params object.
    '''
    
    # Calculate the length of the streams (km)
    params.riverlength = [params.L1]
    for iorder in range(1,params.norder):
        params.riverlength.append(params.riverlength[-1]*params.Rl)

    # Calculate the area of the streams (km2)
    params.riverarea = [params.A1]
    for iorder in range(1,params.norder):
        params.riverarea.append(params.riverarea[-1]*params.Ra)

    # Calculate the number of the streams, first create the list and after that fill the list.
    params.number_of_rivers = []
    for iorder in range(params.norder):
        params.number_of_rivers.append(1.)

    for iorder in range(params.norder-2,-1,-1):
        params.number_of_rivers[iorder] = params.number_of_rivers[iorder+1]*params.Rb

    # Calculate the distribution of N or P load for each stream order as porpotion river*length of the total riverlength
    total_length = []
    for iorder in range(params.norder):
        total_length.append(params.riverlength[iorder] * params.number_of_rivers[iorder])
    params.load_distribution = [] 
    for iorder in range(params.norder):
        params.load_distribution.append(total_length[iorder]/sum(total_length))
     
    # Calculate the flow path probabilities from stream order i to order j as proportion of total length in order j relative to 
    # the total length of all stream orders greater than i. First create a matrix and after that fill the matrix (only the upper part).
    params.flowpath_distribution = []
    for iorder in range(params.norder):
        params.flowpath_distribution.append(params.norder*[0.0])
    for iorder in range(params.norder):
        sumlength = 0.0
        for iorder2 in range(iorder+1,params.norder):
            sumlength += total_length[iorder2]
        for iorder2 in range(iorder+1,params.norder):
            params.flowpath_distribution[iorder][iorder2] = total_length[iorder2]/sumlength

    params.A0 = params.A1/float(params.Ra)


def make_arguments_strahler():
    args_strahler = general_class.General()
    args_strahler.add_item("slope",None)
    args_strahler.add_item("riverlength",None)
    args_strahler.add_item("riverarea",None) 
    args_strahler.add_item("number_of_rivers",None) 
    args_strahler.add_item("load_distribution",None)
    args_strahler.add_item("flowpath_distribution",None) 
    return args_strahler


def define_subgrid_streamorder_riverbasin(params,mask):
    '''
    This function defines the gridded factors for the subgrid_streamorder.
    All output parameters are set in a general class argument object.
    The stream orders 1 to params.norder are defined, 
    where the last order is the first stream that is modelled by PCRGLOBWB.
    '''
    
    # Extract characteristics of the streams from ASCII files.
    # Each characteristic is stored in a dictionnary.
    # Keys are local cell numbers.
    slope = {}
    riverlength = {}
    riverarea = {}
    number_of_rivers = {}
    load_distribution = {}
    flowpath_distribution = {}
    args_strahler = {}
    # Initialize every cell with an empty list
    for icell in range(len(mask)):
        slope[icell] = []
        riverlength[icell] = []
        riverarea[icell] = []
        number_of_rivers[icell] = []
        load_distribution[icell] = []
        flowpath_distribution[icell] = []
        args_strahler[icell] = make_arguments_strahler()
        
    # Loop over all subgrid orders, to get all information.
    for iorder in range(params.norder):

        # Read all input files.
        slope_grid = ascraster.Asciigrid(ascii_file=\
                           os.path.join(params.strahlerdir,"slope_order"+str(iorder+1)+".asc"),\
                           mask=mask,numtype=float)
        riverlength_grid = ascraster.Asciigrid(ascii_file=\
                           os.path.join(params.strahlerdir,"length_order"+str(iorder+1)+".asc"),\
                           mask=mask,numtype=float)
        riverarea_grid = ascraster.Asciigrid(ascii_file=\
                           os.path.join(params.strahlerdir,"basinarea_order"+str(iorder+1)+".asc"),\
                           mask=mask,numtype=float)
        nb_grid = ascraster.Asciigrid(ascii_file=\
                           os.path.join(params.strahlerdir,"nb_order"+str(iorder+1)+".asc"),\
                           mask=mask,numtype=float)
        loaddistrib_grid = ascraster.Asciigrid(ascii_file=\
                           os.path.join(params.strahlerdir,"load_distrib_order"+str(iorder+1)+".asc"),\
                           mask=mask,numtype=float)

        flowdistrib_grids = []
        for j in range(iorder+1,params.norder):
            flowdistrib_grids.append(ascraster.Asciigrid(ascii_file=\
                           os.path.join(params.strahlerdir,"flow_distrib_order"+str(iorder+1)+"_"+str(j+1)+".asc"),\
                           mask=mask,numtype=float))

        # TODO: add option if missing values, e.g. arctic areas
        # Put all info from the grid files into the cell data. So every cell has info for all orders in one list.
        for icell in range(len(mask)):
            slope[icell].append(slope_grid.get_data(icell))
            riverlength[icell].append(riverlength_grid.get_data(icell))
            riverarea[icell].append(riverarea_grid.get_data(icell))
            number_of_rivers[icell].append(nb_grid.get_data(icell))
            load_distribution[icell].append(loaddistrib_grid.get_data(icell))
            flowpath_distribution[icell].append(params.norder*[0.0])
            for j in range(iorder+1,params.norder):
                j2 = j-(iorder+1)
                flowpath_distribution[icell][iorder][j] = flowdistrib_grids[j2].get_data(icell)

    # Set this parameter A0 and put it into the params class.
    params.A0 = params.A1/float(params.Ra)
    
    # Put everythin into one data structure.
    for icell in range(len(mask)):
        args_strahler[icell].set_val("slope",slope[icell])
        args_strahler[icell].set_val("riverlength",riverlength[icell])
        args_strahler[icell].set_val("riverarea",riverarea[icell])
        args_strahler[icell].set_val("number_of_rivers",number_of_rivers[icell])
        args_strahler[icell].set_val("load_distribution",load_distribution[icell])
        args_strahler[icell].set_val("flowpath_distribution",flowpath_distribution[icell])
 
    return args_strahler


