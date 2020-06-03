# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/mouth.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import write_dict
import ascraster
import my_sys
from error import *

import copy

class Mouth:
    '''General functions for the classes NMouth and PMouth'''
    def set_index_grid(self,grid):
        if grid.masked:
            # Calculate index in an unmasked grid
            irow = grid.nrows - int((self.lat-float(grid.yllcorner))/float(grid.cellsize)) - 1
            icol = int((self.lon-float(grid.xllcorner))/float(grid.cellsize))
            ind = irow * int(grid.ncols) + icol
            try:
                self.index_grid = grid.mask.index(ind)
            except ValueError:
                # No grid cell found in the mask with this coordinates.
                print("No grid cell found in the mask with the coordinates "+\
                      str(self.lat) + " " + str(self.lon))
                self.index_grid = -1
        else:
            self.index_grid = grid.get_index_from_coordin(self.lon,self.lat)
    
    def get_index_grid(self,grid):
        if grid.masked:
            return self.index_grid
        else:    
            return self.index_grid_full

    def weigh_avg(self,m,avg_attr,weigh_attr):
        '''
        Calculate the weighed average of two Mouth objects. 
        avg_attr and weigh_attr must be attributes of Mouth.
        '''
        try:
            sum_avg = getattr(self,avg_attr) * getattr(self,weigh_attr) +\
                      getattr(m,avg_attr) * getattr(m,weigh_attr) 
            sum_weigh = getattr(self,weigh_attr) + getattr(m,weigh_attr) 
            setattr(self,avg_attr,sum_avg/sum_weigh)
        except ZeroDivisionError:
            setattr(self,avg_attr,0.0)

    def add_item1(self,name,initvalue = 0.0):
        '''
        Put a parameter in the class with the same name and header.
        '''
        self.add_item(name,header=name,initvalue = initvalue)

    def add_item(self,name,header=None,initvalue = 0.0):
        '''
        Put a parameter in the class with header and initial value.
        Only parameters which have a header, will written to output file in the write function.
        '''
        # Check whether this parameter is in the list.
        if (name in self.names):
            raise MyError("Parameter " + name + " is already in this class.")

        # Put parameter in class (names and headers).
        self.names.append(name)
        self.headers.append(header)

        # Create the attribute name
        setattr(self,name,initvalue)


    def set_item(self,name,val):
        '''
        Set value of parameter name
        '''
        # Check whether this parameter is in the list.
        try:
            qq = getattr(self,name)
        except AttributeError:
            raise MyError("Can not set parameter. Parameter " + name + " is not in this class.")

        # set parameter
        setattr(self,name,val)

    def get_item(self,name):
        '''
        Get value of parameter name
        '''
        # Check whether this parameter is in the list.
        try:
            return getattr(self,name)
        except AttributeError:
            raise MyError("Parameter " + name + " is not in this class.")     
        
    def write_header(self, fp, sep=","):
        '''
        Write the header that belongs to the function write (see below) to the file which has filepointer fp.
        '''
        
        for header in self.headers:
            if (header != None):
                fp.write(header + sep)
        fp.write("\n")
        
    def write(self, fp, sep=","):
        '''
        Write this class to file which has filepointer fp.
        '''
        for item in range(len(self.headers)):
            if (self.headers[item] != None):
                fp.write(str(self.get_item(self.names[item])) + sep)
        fp.write("\n")

    def sum_attr(self,obj,name):
        '''
        Sum the same attribute of two objects in self.
        '''
        self.set_item(name,self.get_item(name) + obj.get_item(name))

 
def read_file(filename,sep,my_object,keyitem="basinid"):
    '''
    Read the file and put everything in a dictionary with 
    int(basinid) as key and a instance of a class (PMouth or NMouth) as value 
    '''
    dict = {}
    # Open file and read all the lines    
    lines = my_sys.my_readfile(filename,my_sys_ldebug=0,fatal=1)
    
    # The header gives the name of the column, which is the same as the item of the class.
    header = lines[0].split(sep)

    # Find key:
    ikey = -1
    for item in range(len(header)):
        if (header[item] == keyitem):
            ikey = item
    print(ikey, header[ikey])

    for line in range(1,len(lines)):
        fields = lines[line].split(sep)
        try:
            key = int(fields[ikey])
        except ValueError:
            key = fields[ikey]
 
        dict[key] = my_object.copy()
        # Put key element in the class
        setattr(dict[key],keyitem,key)
        for item in range(min(len(fields),len(header))):
            if (item != ikey):
                # Put all the values in the Mouth object
                try:
                    dum = int(fields[item])
                    setattr(dict[key],header[item],int(fields[item]))
                except ValueError:
                    try:
                        dum = float(fields[item])
                        setattr(dict[key],header[item],float(fields[item]))
                    except ValueError: 
                        setattr(dict[key],header[item],str(fields[item]))
            
    return dict

def aggregate(dict,name_id_grid,mask,filename,sep,keyname="basinid"):    
    '''
    Aggreggate a dictionary of Mouth objects to another id grid and write the aggregation to output file.
    '''
    
    # Make a copy of the dictionary for the output
    new_dict={}

    # Read id_grid
    id_grid = ascraster.Asciigrid(ascii_file=name_id_grid,mask=mask,numtype=int)

    for key in list(dict.keys()):
        ind = dict[key].index_grid
        if (ind == None):
            # There is no spatial location in isogrid found for this class.
            id = -9999
        else:
            id = int(id_grid.get_data(ind,-9999))
        try:
            new_dict[id].sum(dict[key])
        except KeyError:
            new_dict[id] = dict[key].copy()

    # Make a correction on the name
    for key in list(new_dict.keys()):
        setattr(new_dict[key],keyname,key)
    
    # Aggregate to global level
    for key in list(new_dict.keys()):
        try:
            new_dict["Total"].sum(new_dict[key])
        except KeyError:
            new_dict["Total"] = new_dict[key].copy()

    # Write all information of the river basin to output file:    
    fp = open(filename,"w")
    for key in list(new_dict.keys()):
        new_dict[key].write_header(fp,sep=sep)
        break
    write_dict.write_dict(fp_in=fp,filename=None,dic=new_dict,headerkey="",sep=sep,lkey=0,ldebug=1,fatal=0)
    fp.close()
