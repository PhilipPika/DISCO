# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/make_dict_lake.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import ascraster

class Lake:
    def __init__(self):
       pass

    def write_header(self,fp, sep=" "):
        '''
        Write the header that belongs to the function write (see below) to the file which has filepointer fp.
        '''
        fp.write(str("lakeid") + sep +\
              str("basinid") + sep +\
              str("lendolake") + sep +\
              str("number_of_cells") + sep +\
              str("waterarea") + sep +\
              str("outflow_cell") + sep +\
              str("volume") + sep +\
              str("discharge") + sep +\
              str("depth") + sep +\
              str("hydraulic_load") + sep +\
              str("residence_time")  + sep +\
              str("temperature")  + sep +\
              str("fraction_water")  + sep +\
              str("load_in")  + sep +\
              str("load_out")  + sep +\
              str("retention_load")  + sep +\
              str("retention") + "\n")

    def write(self,fp, sep=" "):
        '''
        Write this class to file which has filepointer fp.
        '''
        fp.write(str(self.id) + sep +\
              str(self.basinid) + sep +\
              str(self.lendolake) + sep +\
              str(len(self.cells)) + sep +\
              str(self.water_area) + sep +\
              str(self.outcell) + sep +\
              str(self.volume) + sep +\
              str(self.discharge) + sep +\
              str(self.depth) + sep +\
              str(self.hydraulic_load) + sep +\
              str(self.residence_time))
        try:
            fp.write(sep + str(self.temperature))
        except AttributeError:
            fp.write(sep +"-9999.")
        try:
            fp.write(sep + str(self.fraction_water))
        except AttributeError:
            fp.write(sep +"-9999.")
        try:
            fp.write(sep + str(self.load_in))
        except AttributeError:
            fp.write(sep +"-9999.")
        try:
            fp.write(sep + str(self.load_out))
        except AttributeError:
            fp.write(sep +"-9999.")
        try:
            inload = float(self.load_in)
            outload = float(self.load_out)
            fp.write(sep + str(inload - outload))
        except AttributeError:
            fp.write(sep +"-9999.")        
        try:
            fp.write(sep + str(self.retention))
        except AttributeError:
            fp.write(sep +"-9999.")
        fp.write("\n")

def calculate(params,mask,discharge_grid):
    '''
    Make a dictionary for the lake properties.
    Take as key the cell number for the dictionary lake_icell_dict.
    Take as key the lakeid for the dictionary lakeid_dict.
    '''
    # Conversion from m to km
    m2tokm2 = 1.0e-6
    m3tokm3 = 1.0e-9
    kmtom   = 1.0e3

    # Define empty dictionaries. 
    lakeid_dict = {}
    lake_icell_dict = {}

    # Read input grids.
    lakeid = ascraster.Asciigrid(ascii_file=params.lakeid,mask=mask,numtype=int)
    endo_lakes = ascraster.Asciigrid(ascii_file=params.endo_lakes,mask=mask,numtype=int)
    outlakeid = ascraster.Asciigrid(ascii_file=params.outlakeid,mask=mask,numtype=int)
    water_area = ascraster.Asciigrid(ascii_file=params.water_area,mask=mask,numtype=float)
    water_storage = ascraster.Asciigrid(ascii_file=params.water_storage,mask=mask,numtype=float)

    for icell in range(lakeid.length):
        id = int(lakeid.get_data(icell,-1))
        if (id > 0):
           # We have found a lake.
           lake_icell_dict[icell] = Lake()
           setattr(lake_icell_dict[icell],"id",id)
           setattr(lake_icell_dict[icell],"cells",[])
           # Check whether we have this lake already found:
           try:
               dum = lakeid_dict[id]
               # Add this list to the dictionary
               lakeid_dict[id].cells.append(icell)
           except KeyError:
               # This one is new.
               lakeid_dict[id] = Lake()
               setattr(lakeid_dict[id],"id",id)
               # Make a list of cell numbers which belong to this lake
               setattr(lakeid_dict[id],"cells",[icell])
               # Set outstream cell of this lake
               setattr(lakeid_dict[id],"outcell",None)
               # Set waterarea of this lake
               setattr(lakeid_dict[id],"water_area",None)
               # Set logical whether this lake has an outstream (endoreic)
               setattr(lakeid_dict[id],"lendolake",None)
               # Set volume of this lake
               setattr(lakeid_dict[id],"volume",None)
               # Set depth of this lake
               setattr(lakeid_dict[id],"depth",None)
               # Set hydraulic load of this lake
               setattr(lakeid_dict[id],"hydraulic_load",None)
               # Set residence time of this lake
               setattr(lakeid_dict[id],"residence_time",None)


           # Check whether cell is outstream point
           val = outlakeid.get_data(icell,-1)
           if (val > 0):
               setattr(lakeid_dict[id],"outcell",icell)            
               setattr(lake_icell_dict[icell],"outcell",1)
               # Volume in km3
               volume = water_storage.get_data(icell,-1)* m3tokm3
               setattr(lakeid_dict[id],"volume",volume)            
               setattr(lake_icell_dict[icell],"volume",volume)
           else:
               setattr(lake_icell_dict[icell],"outcell",0)

           # Check whether cell is endoreic
           val = endo_lakes.get_data(icell,-1)
           if (val > 0):
               setattr(lakeid_dict[id],"lendolake",1)            
               setattr(lake_icell_dict[icell],"lendolake",1)
           else:
               setattr(lakeid_dict[id],"lendolake",0)
               setattr(lake_icell_dict[icell],"lendolake",0)

           # Get total water area of this object
           val = water_area.get_data(icell,-1)
           if (val > 0):
                # Water area in km2
               setattr(lakeid_dict[id],"water_area",val*m2tokm2)            
               setattr(lake_icell_dict[icell],"water_area",val*m2tokm2)
           else:
               setattr(lake_icell_dict[icell],"water_area",0.0)

    # Make for all the lakes a outflow cell
    for key in list(lakeid_dict.keys()):
        icell_outflow =  lakeid_dict[key].outcell
        if (icell_outflow == None):  
            # Make one of the cell of an endoreic lake outflow_cell. Look for the the cell with most water_storage
            icell_outflow = -1
            ldouble = False
            for icell in lakeid_dict[key].cells:
                volume = water_storage.get_data(icell,-1)* m3tokm3
                if (volume > 0.):
                    if (icell_outflow < 0):
                        icell_outflow = icell
                        pot_cells = [icell]
                        vol_cells = [volume]
                    else:
                        if (lakeid_dict[key].lendolake == 1):
                            print("****ERROR****: Endoreic lake ",key," found with more than one cell with a positive water volume.")
                        ldouble = True
                        pot_cells.append(icell)
                        vol_cells.append(volume)
            if (ldouble):
                vol_outflow = vol_cells[0]
                pot_outflow = pot_cells[0]
                for item in range(1,len(pot_cells)):
                    if (vol_cells[item] > vol_outflow):
                        vol_outflow = vol_cells[item]
                        pot_outflow = pot_cells[item]
                icell_outflow = pot_outflow

            if (icell_outflow == -1):
                print("****ERROR****: No outlet found with positive water volume for lake: ",key, " Lake removed.")
                # No water volume found in this lake, so no waterbody, so no lake.
                # Remove from dictionary
                for item in range(len(lakeid_dict[key].cells)):
                    del lake_icell_dict[lakeid_dict[key].cells[item]]
                del lakeid_dict[key]
            else:
                # Set all the attributes of the lakes.
                lakeid_dict[key].outcell = icell_outflow
                lake_icell_dict[icell_outflow].outcell = 1
                volume = water_storage.get_data(icell_outflow,-1)* m3tokm3
                setattr(lakeid_dict[key],"volume",volume)

    # Check whether there are lakes with two or more outlets.
    for icell in lake_icell_dict:
        if (lake_icell_dict[icell].outcell == 1):
            if (lakeid_dict[lake_icell_dict[icell].id].outcell != icell):
                print("****ERROR****: Lake/reservoir ", lake_icell_dict[icell].id," has more than one outflow cell.")
                # Correct this for the time being.
                key = lake_icell_dict[icell].id
                icell_outflow = lakeid_dict[key].outcell
                vol = lakeid_dict[key].volume
                for icell in lakeid_dict[key].cells:
                    vol1 = water_storage.get_data(icell,-1)* m3tokm3
                    if (vol1 > vol):
                        icell_outflow = icell
                        vol = vol1
                # Set the dictionary right
                lakeid_dict[key].outcell = icell_outflow
                lakeid_dict[key].volume = vol
                lake_icell_dict[icell_outflow].outcell = 1
                for icell in lakeid_dict[key].cells:
                    if (icell != icell_outflow):
                        lake_icell_dict[icell].outcell = 0

    # Add the cells attribute also to the outflow points.
    for key in list(lakeid_dict.keys()):
        icell_outflow =  lakeid_dict[key].outcell
        for icell in lakeid_dict[key].cells:
            lake_icell_dict[icell_outflow].cells.append(icell)
        # Add depth and volume to the dictionary
        vol = getattr(lakeid_dict[key],"volume")
        setattr(lake_icell_dict[icell_outflow],"volume",vol)
        swarea = getattr(lakeid_dict[key],"water_area")
        if (swarea == None):
            # There is no lake/reservoir
            print("****ERROR****: No outlet found with positive water area for lake: ",key, " Lake removed.")
            # No water area found in this lake, so no waterbody, so no lake.
            # Remove from dictionary
            for item in range(len(lakeid_dict[key].cells)):
                del lake_icell_dict[lakeid_dict[key].cells[item]]
            del lakeid_dict[key]
            continue                       
        try:
            # Depth in meters
            depth = kmtom * (vol/swarea)
        except ZeroDivisionError:
            depth = 0.0 
        setattr(lake_icell_dict[icell_outflow],"depth",depth)
        setattr(lakeid_dict[key],"depth",depth)
        discharge = discharge_grid.get_data(icell_outflow,0.0)
        setattr(lakeid_dict[key],"discharge",discharge)
        try:
            # Residence time in years    
            residence_time = vol/discharge
        except ZeroDivisionError:
            residence_time = params.max_residence_time_per_lake 
        setattr(lake_icell_dict[icell_outflow],"residence_time",residence_time)
        setattr(lakeid_dict[key],"residence_time",residence_time)
        try:
            # Hydraulic load in m/yr    
            hydraulic_load = lakeid_dict[key].depth/lakeid_dict[key].residence_time
        except ZeroDivisionError:
            hydraulic_load = 0.0 
        setattr(lake_icell_dict[icell_outflow],"hydraulic_load",hydraulic_load)
        setattr(lakeid_dict[key],"hydraulic_load",hydraulic_load)

    print("There are ", len(lakeid_dict), " lakes/reservoirs found in the files.")
    print("There are ", len(lake_icell_dict), " lake/reservoir cells found in the files.") 
    sum1 = 0.0
    for key in lakeid_dict:
        sum1 += lakeid_dict[key].water_area
    print("Water area of lakes and reservoirs is " + str(sum1))

    # Assign river id where lake lies.
    basin = ascraster.Asciigrid(ascii_file=params.basin,mask=mask,numtype=int)
    for key in list(lakeid_dict.keys()):
        icell_outflow =  lakeid_dict[key].outcell
        id = int(basin.get_data(icell_outflow,-9999))
        setattr(lakeid_dict[key],"basinid",id)

    return lakeid_dict,lake_icell_dict
   
