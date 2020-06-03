# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/litho_class.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

from error import *
import ascraster
import os


# Litho PO4(mg P/l) Ea(kJ/mol)     DT50_shallow[y] porosity   Description  
# 0     0.000           0.0        0.0             0.00       Dummy class
# 1     0.000           0.0        0.0             0.00       Major Water bodies (WB)
# 2     0.000           0.0        0.0             0.00       Ice + Glaciers (IG)
# 3     0.142       50000.0        5.0             0.05       Plutonic basic (PB)
# 4     0.018       60000.0        5.0             0.02       Plutonic acid (PA)
# 5     0.085       50000.0        5.0             0.05       Volcanic basic (VB)
# 6     0.008       60000.0        5.0             0.05       Volcanic acid (VA)
# 7     0.020       60000.0        5.0             0.02       Precambrian Basement (PR)
# 8     0.012       60000.0        5.0             0.02       Metamorphic Rocks (MT)
# 9     0.016       60000.0        5.0             0.02       Complex lithology (CL)
# 10    0.015       60000.0        1.0             0.10       Silici-clastic consolidated sedimentary (SS)
# 11    0.018       60000.0        5.0             0.10       Mixed consolidated sedimentary (SM)
# 12    0.073           0.0        5.0             0.10       Carbonated consolidated sedimentary (SC)
# 13    0.000           0.0        5.0             0.20       Evaporites (EP)
# 14    0.017       60000.0        5.0             0.30       Semi- to unconsolidated sedimentary (SU)
# 15    0.017       60000.0        2.0             0.15       Alluvial deposits (AD)
# 16    0.010       60000.0        5.0             0.20       Loess (LO)
# 17    0.010       60000.0        5.0             0.30       Dunes and shifting sand (DS)


class Litho:
    '''Derived parameter values of the lithology map'''
    def __init__(self,backgroundPO4=0.0,dt50_shallow=0.,dt50_deep=0.,porosity=0.,permeable=0,bSi_weathering=0.0,Ea=0.0):
        self.backgroundPO4 = float(backgroundPO4)
        self.dt50_shallow = float(dt50_shallow)
        self.dt50_deep = float(dt50_deep)
        self.porosity = float(porosity)
        self.bSi_weathering = float(bSi_weathering)
        self.permeable = float(permeable)
        # Activation energy
        self.Ea = float(Ea)

    def make_value_list(self,params,mask,list_attribute_names):
        '''
        The litho grid is given in classes. For each cell we know the fraction of that class in the cell.
        Grid0 is the nodata information. When grid0 is 1.0 then the cell is 100% nodata.
        A list with length of grids (normally length of mask) is returned with on each entry the
        value of the attribute given in the list_attribute_names (and also the order as given in this list).
        '''
        # Compare first grid with another grid of the run
        equal = ascraster.compare_grids(params.pgrass,os.path.join(params.litho_dir,"nodatgrid.asc"))
        if (not equal):
            raise MyError("Litho maps don't have same header as other input files.")
        
        grid0 = ascraster.Asciigrid(ascii_file=os.path.join(params.litho_dir,"nodatgrid.asc"),\
                                       mask=mask,numtype=float)

        # Check whether there are fractions given in grid0
        max_value = max(grid0.values)
        if (max_value > 1.0):
            raise MyError("Litho maps must be a fraction. Found values greater than 1 (found: " + str(max_value) + ").")

    
        # Make a default output list:
        out = []
        for i in range(grid0.length):
            out.append(len(list_attribute_names) * [0.0])
        
        # Loop over all the classes of litho map. Start with class 3 because first two are major water bodies and ice:
        for iclass in range(3,18):
            # Get the litho_class values for this list_attributes_names
            values = []
            for item in range(len(list_attribute_names)):
                values.append(float(getattr(litho_lib[iclass],list_attribute_names[item])))

            # Check header of input file
            equal = ascraster.compare_grids(params.pgrass,os.path.join(params.litho_dir,"grid" + str(iclass) + ".asc"))
            if (not equal):
                raise MyError("Litho maps don't have same header as other input files.")                

            grid = ascraster.Asciigrid(ascii_file=os.path.join(params.litho_dir,"grid" + str(iclass) + ".asc"),\
                                       mask=mask,numtype=float)
           
            for icell in range(0,grid0.length):
                grid0_cell = grid0.get_data(icell,0.0)
                if (grid0_cell < 1.0):
                    # There is something to do.
                    grid_cell = grid.get_data(icell,0.0)

                    if (grid_cell > 0.0):
                        # There is something of this class in this cell:                       
                        for item in range(len(list_attribute_names)):
                            out[icell][item] += grid_cell * values[item]

        correc = grid0.length * [0.0]
        
        # Correction loop to correct for the nodata of the cell
        for iclass in range(1,3):
            # Get the litho_class values for this list_attributes_names

            # Check header of input file
            equal = ascraster.compare_grids(params.pgrass,os.path.join(params.litho_dir,"grid" + str(iclass) + ".asc"))
            if (not equal):
                raise MyError("Litho maps don't have same header as other input files.")                

            grid = ascraster.Asciigrid(ascii_file=os.path.join(params.litho_dir,"grid" + str(iclass) + ".asc"),\
                                       mask=mask,numtype=float)
           
            for icell in range(0,grid0.length):
                grid0_cell = grid0.get_data(icell,0.0)
                if (grid0_cell < 1.0):
                    # There is something to do.
                    if (iclass == 1):
                        correc[icell] += grid0_cell

                    grid_cell = grid.get_data(icell,0.0)
                    if (grid_cell > 0.0):
                        correc[icell] += grid_cell

        # Make correction for the situation where there is water bodies or ice or a part of the cell is nodata.
        for icell in range(len(correc)):
            if (correc[icell] > 0.0):
                # There was ice or water bodies or part of the cell is nodata

                if (correc[icell] < 1.0):
                    for item in range(len(list_attribute_names)):
                        out[icell][item] /= (1.0 - correc[icell])
        return out
        
# Si: Moosdorf et al (2010) Calibration using the Durr et al (2005).

# Make a database with all characteristics into one library
# DT50_deep is zero for all litho classes.
litho_lib = []

# Litho PO4(mg P/l) Ea(kJ/mol)     DT50_shallow[y] porosity   Description
# 0     0.0000          0.0        0.0             0.00       Dummy class
# Here bSi_weathering is the value of the intercept.
litho_lib.append(Litho(backgroundPO4=0.000,Ea=0.0,dt50_shallow=0.0,dt50_deep=0.0,porosity=0.00,permeable=0,bSi_weathering=0.0216))

# 1     0.0000          0.0        0.0             0.00       Major Water bodies (WB)
litho_lib.append(Litho(backgroundPO4=0.000,Ea=0.0,dt50_shallow=0.0,dt50_deep=0.0,porosity=0.00,permeable=0,bSi_weathering=0.0))

# 2     0.0000          0.0        0.0             0.00       Ice + Glaciers (IG)
litho_lib.append(Litho(backgroundPO4=0.000,Ea=0.0,dt50_shallow=0.0,dt50_deep=0.0,porosity=0.00,permeable=0,bSi_weathering=0.0))

# 3     0.142       50000.0        5.0             0.05       Plutonic basic (PB)
# Si: Comes from Jansen (= Moosdorf) et al (2010)
litho_lib.append(Litho(backgroundPO4=0.142,Ea=50000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.05,permeable=0,bSi_weathering=1.01))

# 4     0.018       60000.0        5.0             0.02       Plutonic acid (PA)
litho_lib.append(Litho(backgroundPO4=0.018,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.02,permeable=0,bSi_weathering=0.9845))

# 5     0.085       50000.0        5.0             0.05       Volcanic basic (VB)
litho_lib.append(Litho(backgroundPO4=0.085,Ea=50000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.05,permeable=0,bSi_weathering=0.9501))

# 6     0.008       60000.0        5.0             0.05       Volcanic acid (VA)
# Si: Comes from Jansen (= Moosdorf) et al (2010)
litho_lib.append(Litho(backgroundPO4=0.008,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.05,permeable=0,bSi_weathering=0.8191))

# 7     0.020       60000.0        5.0             0.02       Precambrian Basement (PR)
litho_lib.append(Litho(backgroundPO4=0.020,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.02,permeable=0,bSi_weathering=0.8593))

# 8     0.012       60000.0        5.0             0.02       Metamorphic Rocks (MT)
litho_lib.append(Litho(backgroundPO4=0.012,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.02,permeable=0,bSi_weathering=0.8771))

# 9     0.016       60000.0        5.0             0.02       Complex lithology (CL)
# Si: Same value as SS
litho_lib.append(Litho(backgroundPO4=0.016,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.02,permeable=0,bSi_weathering=0.8389))

# 10    0.015       60000.0        1.0             0.10       Silici-clastic consolidated sedimentary (SS)
litho_lib.append(Litho(backgroundPO4=0.015,Ea=60000.0,dt50_shallow=1.0,dt50_deep=0.0,porosity=0.10,permeable=0,bSi_weathering=0.8389))

# 11    0.018       60000.0        5.0             0.10       Mixed consolidated sedimentary (SM)
# Si: Comes from Jansen (= Moosdorf) et al (2010)
litho_lib.append(Litho(backgroundPO4=0.018,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.10,permeable=0,bSi_weathering=0.7531))

# 12    0.073          0.0        5.0             0.10       Carbonated consolidated sedimentary (SC)
litho_lib.append(Litho(backgroundPO4=0.073,Ea=    0.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.10,permeable=0,bSi_weathering=0.7827))

# 13    0.0000          0.0        5.0             0.20       Evaporites (EP)
litho_lib.append(Litho(backgroundPO4=0.0000,Ea=    0.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.20,permeable=0,bSi_weathering=0.0))

# 14    0.017       60000.0        5.0             0.30       Semi- to unconsolidated sedimentary (SU)
litho_lib.append(Litho(backgroundPO4=0.017,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.30,permeable=1,bSi_weathering=0.7961))

# 15    0.017       60000.0        2.0             0.15       Alluvial deposits (AD)
# Si: Same value as SU
litho_lib.append(Litho(backgroundPO4=0.017,Ea=60000.0,dt50_shallow=2.0,dt50_deep=0.0,porosity=0.15,permeable=1,bSi_weathering=0.7961))

# 16    0.010      60000.0        5.0             0.20       Loess (LO)
# Si: Same value as SU
litho_lib.append(Litho(backgroundPO4=0.010,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.20,permeable=1,bSi_weathering=0.7961))

# 17    0.010      60000.0        5.0             0.30       Dunes and shifting sand (DS)
# Si: Same value as SU
litho_lib.append(Litho(backgroundPO4=0.010,Ea=60000.0,dt50_shallow=5.0,dt50_deep=0.0,porosity=0.30,permeable=1,bSi_weathering=0.7961))


def calculate_water_litho(params,mask):
    '''
    Calculate the total discharge for each litho class. 
    This is needed to fill in the table of litho_lib with information from Jens Hartmann
    '''
    pnet = ascraster.Asciigrid(ascii_file=params.pnet,mask=mask,numtype=float)
    if params.ldebug: print(params.pnet  + " has been read.")  

    landarea = ascraster.Asciigrid(ascii_file=params.landarea,mask=mask,numtype=float)
    if params.ldebug: print(params.landarea  + " has been read.")    
    
    # Calculate accumulated water flux of each grid cell to the mouth of the river basin. Make discharge in km3
    # Convert pnet [mm/yr] to km3 by multiplying with landarea and 10e-6 (conversion from mm to km)
    water = ascraster.duplicategrid(landarea)    
    for icell in range(landarea.length):
        wt = pnet.get_data(icell,0.0) * landarea.get_data(icell,0.0) * 1.0e-6
        water.set_data(icell,wt)

    # Compare first grid with another grid of the run
    equal = ascraster.compare_grids(params.pgrass,os.path.join(params.litho_dir,"nodatgrid.asc"))
    if (not equal):
        raise MyError("Litho maps don't have same header as other input files.")
        
    grid0 = ascraster.Asciigrid(ascii_file=os.path.join(params.litho_dir,"nodatgrid.asc"),\
                                mask=mask,numtype=float)

    # Check whether there are fractions given in grid0
    max_value = max(grid0.values)
    if (max_value > 1.0):
        raise MyError("Litho maps must be a fraction. Found values greater than 1 (found: " + str(max_value) + ").")
  
    # Make a default output list:
    out = len(litho_lib) * [0.0]
    area = len(litho_lib) * [0.0]

    # Loop over all the classes of litho map. Start with class 3 because first two are major water bodies and ice:
    for iclass in range(1,18):
        print(iclass)
        # Check header of input file
        equal = ascraster.compare_grids(params.pgrass,os.path.join(params.litho_dir,"grid" + str(iclass) + ".asc"))
        if (not equal):
            raise MyError("Litho maps don't have same header as other input files.")                

        grid = ascraster.Asciigrid(ascii_file=os.path.join(params.litho_dir,"grid" + str(iclass) + ".asc"),\
                                   mask=mask,numtype=float)
           
        print(os.path.join(params.litho_dir,"grid" + str(iclass) + ".asc") + " read.")
        for icell in range(0,grid0.length):
            grid0_cell = grid0.get_data(icell,0.0)
            if (grid0_cell < 1.0):
                # There is something to do.
                grid_cell = grid.get_data(icell,0.0)

                if (grid_cell > 0.0):
                    # There is something of this class in this cell:                       
                    out[iclass] += grid_cell * water.get_data(icell,0.0)
                    area[iclass] += grid_cell * landarea.get_data(icell,0.0)
            
            if (iclass == 1):
                out[0]  += grid0_cell * water.get_data(icell,0.0)
                area[0] += grid0_cell * landarea.get_data(icell,0.0)

    return out,area
