# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/ascraster.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
'''
ascraster.py

A Python library containing functions to read and write ASCII raster files.

Initialy created on 26 mei 2009
@author: warrinka
'''
import os
import copy
# For reading compressed files
import itertools
import math

try:
    from PIL import Image
    import ascbitmap
    PILavailable = True
except ImportError:
    PILavailable = False

import pickle as pickle
import gzip
    

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class ASCIIGridError(Error):
    '''Is raised when a ASCII grid file is not valid.
    
    file: file object of the ASCII grid
    message: message with explanation on failure
    '''
    def __init__(self, message, file = None):
        self.message = message
        if file != None:
            self.filename= file.name
            self.__hasfile = True
            # Close the opened file
            file.close()
        else:
            self.__hasfile = False
        
    def __str__(self):
        if self.__hasfile:
            return repr((self.filename, self.message))
        else:
            return repr(self.message)

class Asciigrid(object):
    '''
    Asciigrid reads, manipulates and creates ASCII grids.
    '''
    def __init__(self, ncols = None, nrows = None, xllcorner = None, yllcorner = None, cellsize = None, nodata_value = None,\
                       ascii_file = None, pickle_file = None, mask = None, numtype = None):
        # Set masked property. Will be False until the values are masked
        # masked property is used to make sure a set of values will not be masked twice
        self.masked = False
        self.mask = mask
        # Check value of numtype
        if not (numtype in (float, int, None)):
            if ascii_file: print("ERROR for file: %s" % ascii_file)
            raise Exception("Type of numbers can be float or int or None. Wrong type specified.")
        self.numtype = numtype
        
        # If a ascii_file is given get properties and values from this ascii_file
        if ascii_file:
            # Set header properties
            headerprops = read_header(ascii_file)
            self.ncols = int(headerprops['ncols'])
            self.nrows = int(headerprops['nrows'])
            try:
                self.xllcorner = self._str2num(headerprops['xllcorner'],True)
                self.yllcorner = self._str2num(headerprops['yllcorner'],True)
                self.cellsize = self._str2num(headerprops['cellsize'],True)
                if (headerprops['nodata_value'] != None):
                    self.__nodata_value = self._str2num(headerprops['nodata_value'],True)
                else:
                    self.__nodata_value = None
            except ValueError as vals:
                raise ASCIIGridError(vals)
            if self.__nodata_value == None:
                self.__nodata = False
            else:
                self.__nodata = True
            # Set values property
            self.length = self.ncols * self.nrows
            self.values = self._read_values(ascii_file)
        elif pickle_file:
            self = read_pickle_file(pickle_file)
            
        # If no ascii_file is given write the properties from the given arguments
        else:
            self.ncols = int(ncols)
            self.nrows = int(nrows)
            self.xllcorner = xllcorner
            self.yllcorner = yllcorner
            self.cellsize = cellsize
            self.__nodata_value = nodata_value
            self.length = self.ncols * self.nrows
            # set the values as a list with None values
            self.values = [nodata_value for x in range(self.length)]        
            
            if self.__nodata_value == None:
                self.__nodata = False
            else:
                self.__nodata = True
       
        if self.mask != None:
            self.apply_mask(self.mask)

   # def __del__(self):
   #     print "Delete function of Asciigrid is called."
   #     del self.values
   #     if self.masked: self.mask=None
   #     gc.collect()
    
    def _get_nodata_value(self):
        '''
        Gets value for nodata
        '''
        return self.__nodata_value
    
    def _set_nodata_value(self, newval):
        '''
        Sets value for nodata. If nodata_value is set to a value other than None 
        the previous nodata cells are also changed into the new value.
        If nodata is False the nodata_value is set but the cells are not changed.
        @param newval: New no data value
        @type newval: INTEGER/FLOAT
        '''
        if newval == None:
            self.nodata = False
            self.__nodata_value = newval
        elif self.nodata == True:
            oldval = self.__nodata_value
            # give the new nodata value to the affected cells
            for i, value in enumerate(self.values):
                if value == oldval:
                    self.values[i] = newval
            self.__nodata_value = newval
        else:
            oldval = None
            # give the new nodata value to the affected cells
            for i, value in enumerate(self.values):
                if value == oldval:
                    self.values[i] = newval
            self.__nodata_value = newval
            self.nodata = True
        
          
    def _get_nodata(self):
        '''
        Returns boolean, whether a nodata value is present.
        '''
        return self.__nodata
    
    def _set_nodata(self, status):
        '''
        Sets the status for nodata. If nodata_value is not provided (is None)
        and status is changed to True an error is raised.
        @param status: New no data value
        @type status: BOOLEAN
        '''
        # If new status is True check if nodata_value not is None
        if status and self.__nodata_value is None:
            raise Exception("no data status can not be set to True if no no data value is provided.")
        else:
            self.__nodata = bool(status)
            
    # nodata must never be changed without changing the cell data with no data
    nodata_value = property(_get_nodata_value, _set_nodata_value)
    nodata = property(_get_nodata, _set_nodata)
    
    def _read_values(self, filename):
        '''Read values from a ASCIIgrid. Returns a list with values.'''
        # Open the Ascii Raster file
        if (compressed(filename)):
            f = gzip.open(filename, 'r')
        else:
            f = open(filename, 'r')
        # skip the first 5 lines header info. 6th line can contain header info but also values.
        if self.nodata == False:
            number_of_header_rows = 5
        else:
            number_of_header_rows = 6
        for _i in range(number_of_header_rows):
            line = f.readline()
        # Set the initial value for the return values.
        values = []
        # Start processing the lines in the asciifile
        try:
            lines = f.readlines(100000000)
            while lines:
                for line in lines:
                    if (self.numtype == float):
                        values.extend(list(map(float, line.split())))
                    elif (self.numtype == int):
                        values += list(map(int, line.split()))
                    else:
                        values += list(map(self._str2num, line.split()))
                lines = f.readlines(100000000)
            
            # check if the returned value list contains the number of values described in the header
            if len(values) != self.length: 
                raise ASCIIGridError("Incorrect number of values. %s found; %s expected." % (len(values), self.length), f)
            f.close()
        except ValueError:
            # Try compress file first before giving an error
            f.close()
            values = self._read_values_compress(filename)
        except:
            raise ASCIIGridError("Error in reading file %s ." % filename)

        return values


    def _read_values_compress(self, filename):
        '''Read values from a ASCIIgrid. Returns a list with values.'''
        # Open the Ascii Raster file
        if (compressed(filename)):
            f = gzip.open(filename, 'r')
        else:
            f = open(filename, 'r')
        # skip the first 5 lines header info. 6th line can contain header info but also values.
        if self.nodata == False:
            number_of_header_rows = 5
        else:
            number_of_header_rows = 6
        for _i in range(number_of_header_rows):
            line = f.readline()
        # Set the initial value for the return values.
        values = []
        # Start processing the lines in the asciifile
        lines = f.readlines(100000000)
        while lines:
            for line in lines:
                try:
                    if (self.numtype == float or self.numtype == int):
                        values.extend(list(map(self._str2num_compress, line.split())))
                    else:
                        values += list(map(self._str2num_compress, line.split()))
                except ValueError:
                    raise ASCIIGridError("Incorrect value found. Expected a numerical value(s): \n%s" % (line), f)
            lines = f.readlines(100000000)

        # Make one long list of the values
        values = list(itertools.chain(*values))

        # check if the returned value list contains the number of values described in the header     
        if len(values) != self.length: 
            raise ASCIIGridError("Incorrect number of values. %s found; %s expected." % (len(values), self.length), f)
        f.close()
        return values
    
    def _str2num(self, s,special=None):
        '''
        checks if a string is a float or integer. Returns a float or integer.
        @param s: string with a numerical value
        @type s: STRING
        '''
        try:
            if special == True:
                try:
                    return int(s)
                except ValueError:
                    return float(s)
            elif (self.numtype == float):
                return float(s)
            elif (self.numtype == int):
               try:
                   return int(s)
               except ValueError:
                   return float(s)

            try:
                return int(s)
            except ValueError:
                return float(s)
        except ValueError:
            raise ValueError("Value is not numerical: %s" % s)

    def _str2num_compress(self, s):
        '''
        checks if a string is a float or integer. Returns a float or integer. Handles strings of num*value (compressed format)
        @param s: string with a numerical value
        @type s: STRING
        '''
        # Try integer
        try:
            return [int(s)]
        except ValueError:
            # Try float
            try:
                return [float(s)]
            except ValueError:
                # Try compressed value
                slist = s.split('*')
                try:
                    times = int(slist[0])
                    val = float(slist[1])
                    return times *  [val]
                except:
                    if (len(slist) != 2):
                        raise ASCIIGridError("Value1 is not numerical: %s" % s)
                    else:
                        try:
                            times = int(slist[0])
                        except ValueError:
                            raise ASCIIGridError("Number of values must be integer. Found value: %s in field %s" % (slist[0],s))
                        try:
                            val = float(slist[1])
                        except ValueError:
                            raise ASCIIGridError("Value is not numerical. Found value: %s in field %s" % (slist[1],s))

                    # For all other situations
                    raise ValueError("Value is not numerical: %s" % s)
        

    def get_data(self, ind,val=None):
        '''
        Get data value for a given index position. When the index position is NoData None is returned.
        When val is not None, val will be returned in stead of None (in case of NoData) 
        Value returned as FLOAT or INTEGER.
        @param ind: Index for the asked value
        @type ind: INTEGER
        @param val: Default value in stead of None in case of Nodata
        @type val: FLOAT
        
        '''
        # If the given index falls outside the given domain, raise ASCIIGridError
        if ind >= self.length or ind < 0: 
            raise ASCIIGridError("Index position %s falls outside bounds." % ind)

        returnValue = self.values[ind]

        # Check if the returned value is the no data value, if so return None
        if returnValue == self.__nodata_value:
            if val == None:
                returnValue = None
            else:
               returnValue = val
               
        return returnValue
        
    def set_data(self, ind, value):
        '''
        Get data value for a given index position.
        @param ind: Index for the asked value
        @type ind: INTEGER
        @param value: New value as float or integer
        @type value: FLOAT/INTEGER
        '''
        if ind >= self.length or ind < 0: 
            raise ASCIIGridError("Index position %s falls outside bounds." % ind)

        self.values[ind] = value
    
    def apply_mask(self, mask):
        '''
        Apply a mask on a instance which is not yet masked. The mask object is a list containing indexes to the cells. This list is
        The mask object is stored in the instance of Asciigrid as a pointer. Therefore the mask object should not be altered.
        @param mask: A list with indexes which will not be masked
        @type mask: LIST
        '''
        if not self.masked:
            # Do not deepcopy mask to self.mask. This saves memory. Downside is that the mask object given should not be altered.
            self.mask = mask
            self.values = [self.values[int(x)] for x in self.mask]
            self.masked = True
            self.length = len(self.values)
        else: 
            raise ASCIIGridError("The raster is already masked.")
            
    def add_values(self, new_values):
        '''
        This function adds new values to the ASCII raster instance.
        The given values must have a length equal to self.length.
        @param new_values: A list with values.
        @type new_values: LIST
        '''

        if len(new_values) == self.length:
            # self.values = copy.deepcopy(new_values)
            self.values = new_values[:]
        else:
            raise ASCIIGridError("Values not valid. A length of %s was expected" % self.length)
        
    def compare(self, filename):
        '''Compares the grid with another grid.
        Returns True if ncols, nrows, xllcorner, yllcorner and cellsize are equal.'''
        equal = True
        # Get props from given ASCIIgrid
        file_headerprops = read_header(filename)
        # Check self props against file props
        if self.ncols != int(file_headerprops['ncols']): equal = False
        if self.nrows != int(file_headerprops['nrows']): equal = False
        if self.xllcorner != self._str2num(file_headerprops['xllcorner'],True): equal = False
        if self.yllcorner != self._str2num(file_headerprops['yllcorner'],True): equal = False
        if self.cellsize != self._str2num(file_headerprops['cellsize'],True): equal = False
        return equal

    def write_ascii_file_slow(self, filename,output_nodata_value=None,compress=False):
        '''Writes a ASCIIgrid to the given location on the file system.
           This routine uses less memory because every line is written to file.
           When memory is not a problem, use write_ascii_file which is much faster.
        '''
        print("Start slow method to write grid to file. ")
        if (compress):
            f = gzip.open(filename+".gz", 'w')
        else:
            f = open(filename, 'w')

        f.write("ncols %s\n" % self.ncols)
        f.write("nrows %s\n" % self.nrows)
        f.write("xllcorner %s\n" % self.xllcorner)
        f.write("yllcorner %s\n" % self.yllcorner)
        f.write("cellsize %s\n" % self.cellsize)

        if (self.masked):
            # A masked grid must have a nodata value specified.
            if (self.nodata_value == None and output_nodata_value == None):
                raise ASCIIGridError(("There is a masked grid without a nodata_value."))
            if (self.nodata_value == None):
                self.nodata_value = output_nodata_value

        nodata = self.nodata
        if (self.nodata_value == None):
            nodata_value = str(output_nodata_value)
            if (nodata_value != str(None)):
                nodata = True
                if (not self.masked):
                    # we have to fill all the cells.
                    self.nodata_value = output_nodata_value
        else:
            nodata_value = str(self.nodata_value)
            nodata = True

        if (nodata):
            f.write("NODATA_value %s\n" % nodata_value)

        # check if a mask is used
        if (self.masked):
            # Use a counter for the mask
            imask = 0
            # iterate over the rows
            for c in range(0, self.ncols*self.nrows, self.ncols):
                line_list = [nodata_value for x in range(self.ncols)]        
                # line_list = list(nodata_value for _i in range(self.ncols))
                for _j in range(imask,self.length):
                    if self.mask[_j] < c+self.ncols:
                        line_list[self.mask[_j]-c] = self.values[_j]
                        imask += 1
                    else:
                        # Go to next row
                        break
                line_list1 = [str(x) for x in line_list]        
                f.write("%s\n" % ' '.join(line_list1))        

            # Check whether every masked item is written
            if imask != self.length:
                raise ASCIIGridError("There is something wrong with writing a mask grid. %s used; %s expected." % (imask, self.length),f)
        else:
            for c in range(0, self.length, self.ncols):        
                f.write("%s\n" % ' '.join(map(str, self.values[c:c+self.ncols])))
        f.close()
        
    def row_column(self,list):
        '''
        Returns a pair of rownumber and columnnumber for a list of index (only unmasked!)
        '''
        return [[int(point)//int(self.ncols), int(point) - int(point)//int(self.ncols)*self.ncols] for point in list]

    def write_ascii_file(self, filename,output_nodata_value=None,compress=False):
        '''Writes a ASCIIgrid to the given location on the file system.
           One write to file of the all the numbers. More memory is needed here.
           If memory is a problem, then use the slow version of this routine (write_ascii_file_slow)
           output_nodata_value is used when the self.nodata_value is None and nodata occurs in self.
        '''
        try:
            if (compress):
                f = gzip.open(filename+".gz", 'w')
            else:
                f = open(filename, 'w')
  
            f.write("ncols %s\n" % self.ncols)
            f.write("nrows %s\n" % self.nrows)
            f.write("xllcorner %s\n" % self.xllcorner)
            f.write("yllcorner %s\n" % self.yllcorner)
            f.write("cellsize %s\n" % self.cellsize)

            if self.masked:
                # A masked grid must have a nodata value specified.
                if (self.nodata_value == None and output_nodata_value == None):
                    raise ASCIIGridError(("There is a masked grid without a nodata_value."))
                if (self.nodata_value == None):
                    self.nodata_value = output_nodata_value

            nodata = self.nodata
            if (self.nodata_value == None):
                nodata_value = str(output_nodata_value)
                if (nodata_value != str(None)):
                    nodata = True
                    if (not self.masked):
                        # we have to fill all the cells.
                        self.nodata_value = output_nodata_value
            else:
                nodata_value = str(self.nodata_value)
                nodata = True

            if nodata:
                f.write("NODATA_value %s\n" % nodata_value)

            # check if a mask is used
            if self.masked:
                line_list = []
                # Make all cells no data
                for c in range(0, self.ncols*self.nrows, self.ncols):
                    line_list.append([nodata_value for x in range(self.ncols)])
                 
                # Mask row and column numbers of the mask entries
                _rc_number = self.row_column(self.mask)
              
                # Fill in the data values of the mask
                for _j in range(0,self.length):
                    line_list[_rc_number[_j][0]][_rc_number[_j][1]] = str(self.values[_j])

                # Join all the strings on one line
                line_list = ["%s\n" % ' '.join(line_list[x]) for x in range(len(line_list))]         
                # Join all the lines
                line = "".join(line_list)
                # Write everything to file
                f.write(line)
            else: 
                for c in range(0, self.length, self.ncols):
                    f.write("%s\n" % ' '.join(map(str, self.values[c:c+self.ncols])))
            f.close()
        except MemoryError:
            print("Not enough memory to write grid to file. Trying slow method.")
            # Clean up some memory and try the slow method.
            if (line_list != None): del line_list
            if (f != None):
                # Close the opened file
                f.close()
            import sys
            sys.exc_clear()
            sys.exc_info()[2] = sys.last_traceback = None
            self.write_ascii_file_slow(filename,\
                  output_nodata_value=output_nodata_value,compress=compress)                       
        
    def save_as_bitmap(self, filename, minV=None, maxV=None,mode=None,color=None, number_of_classes = None, classtype=None,\
                       nodatacolor=None):
        '''
        argument mode can be "L" or "RGB"
        argument color works for mode RGB. Color = (1,1,1) means red, green and blue. If just one or two colors, then make a zero in the tuple.
        '''
        if PILavailable:
            if self.masked:
                list1 = list(self.__nodata_value for _i in range(self.nrows*self.ncols))
                for c in range(self.length):
                    list1[self.mask[c]]=self.values[c]
            else:
                list1 = self.values
            im = ascbitmap.Asciibitmap(mode=mode)
            im.createbitmap(self.nrows, self.ncols, list1, filename, self.__nodata_value, minV=minV, maxV=maxV,color=color,\
                                number_of_classes = number_of_classes, classtype=classtype, nodatacolor=nodatacolor)
        else:
            print("No picture saved.\nPlease install the Python Imaging Library (PIL)")

        
    def row_column(self,list):
        '''
        Returns a pair of rownumber and columnnumber for a list of index (only unmasked!)
        '''
        return [[int(point)//int(self.ncols), int(point) - (int(point)//int(self.ncols))*self.ncols] for point in list]

    def get_coordin_from_row_col(self,irow,icol):
        '''
        Returns mid point of cell given by irow and icol
        '''
        if self.masked:
            raise ASCIIGridError("GET_COORDIN_FROM_ROW_COL: Function can not be used for masked grids")
        x = float(self.xllcorner) + (icol+0.5)*float(self.cellsize)
        y = float(self.yllcorner) + (self.nrows-irow-0.5)*float(self.cellsize)
        return x,y

    def get_coordin_from_index(self,ind):
        '''
        Returns mid point of cell given by irow and icol
        '''
        if self.masked:
            raise ASCIIGridError("GET_COORDIN_FROM_INDEX: Function can not be used for masked grids")        
        irow,icol = self.get_row_col_from_index(ind)
        return self.get_coordin_from_row_col(irow,icol)

    def get_row_col_from_index(self,ind):
        '''
        Returns irow and icol of an index
        '''
        if self.masked:
            raise ASCIIGridError("GET_ROW_COL_FROM_INDEX: Function can not be used for masked grids")
        irow = int(ind)//int(self.ncols)
        icol = int(ind) - irow * self.ncols
        return irow,icol

    def get_row_col_from_coordin(self,x,y):
        '''
        Returns irow and icol of an coordinate
        '''
        if self.masked:
            raise ASCIIGridError("GET_ROW_COL_FROM_COORDIN: Function can not be used for masked grids")
        ind = self.get_index_from_coordin(x,y)
        return self.get_row_col_from_index(ind)

    def get_index_from_coordin(self,x,y):
        '''
        Get index of x- and y-coordinate
        '''      
        if self.masked:
            raise ASCIIGridError("GET_INDEX_FROM_COORDIN: Function can not be used for masked grids")
        irow = self.nrows - int((y-float(self.yllcorner))/float(self.cellsize)) - 1
        #irow = int((float(self.yllcorner) + int(self.nrows)*float(self.cellsize) - y)/float(self.cellsize))
        icol = int((x-float(self.xllcorner))/float(self.cellsize))

        if (icol >= self.ncols or icol < 0):
            raise ASCIIGridError("GET_INDEX_FROM_COORDIN: Y position (%s,%s) falls outside bounds." % (x,y))
        if (irow >= self.nrows or irow < 0):
            raise ASCIIGridError("GET_INDEX_FROM_COORDIN: X position (%s,%s) falls outside bounds." % (x,y))

        return self.get_index_from_row_col(irow,icol)

    def get_index_from_row_col(self,irow,icol):
        '''
        Get index of combination of row and col numbers
        '''
        if self.masked:
            raise ASCIIGridError("GET_INDEX_FROM_ROW_COL: Function can not be used for masked grids")
        ind = irow * int(self.ncols) + icol
        if ind >= self.length or ind < 0: 
            raise ASCIIGridError("GET_INDEX_FROM_ROW_COL: Row,Col position (%s,%s) falls outside bounds." % (irow,icol))
        return ind
    
       
    def get_absolute_length_from_ind(self, ind):
        '''
        Returns latitudinal and longitudinal lengths of gridcell in km 
        '''
        lat = self.get_coordin_from_index(ind)[1]
        eq_radius = 6375.1370
        pol_radius = 6356.7523
        amp = eq_radius-pol_radius
        OP = eq_radius-amp*math.sin((lat/180)*math.pi)
        ZP = OP*math.cos((lat/180)*math.pi)
        lat_peri =  2.*math.pi*ZP
        lon_peri = 2.*math.pi*pol_radius
        lat_length = self.cellsize*lat_peri/360.
        lon_length = lon_peri*(float(self.cellsize)/360)
        return lat_length,  lon_length

    def get_absolute_size_from_ind(self, ind):
        '''
        Returns surface area of gridcell in km2
        '''
        latlen, lonlen = self.get_absolute_length_from_ind(ind)
        return latlen*lonlen
   
    def _multi(self,x,y): 
        return x*y

    def _power(self,x,y):
        try:
            return math.pow(x,y)
        except ValueError:
            return None

    def _substraction(self,x,y): 
        return x-y

    def _summation(self,x,y): 
        return x+y

    def _check_min(self,x,val_min):
        if (x < val_min): 
            return val_min
        else:
            return x

    def _check_max(self,x,val_max):
        if (x > val_max): 
            return val_max
        else:
            return x
     
    def _check_range(self,x,val_min,val_max):
        x1 = self._check_min(x,val_min)
        return self._check_max(x1,val_max)

    def _check_range_values(self,minimum,maximum): 
        # Check the range of the result
        if (minimum != None or maximum != None):
            if (minimum != None and maximum != None):
                self.values[:] = [self._check_range(x,minimum,maximum) for x in self.values]
            elif (minimum != None):
                self.values[:] = [self._check_min(x,minimum) for x in self.values]
            else:
                self.values[:] = [self._check_max(x,maximum) for x in self.values]

    def _check_range_value(self,x,minimum,maximum): 
        # Check the range of the result
        if (minimum != None or maximum != None):
            if (minimum != None and maximum != None):
                return self._check_range(x,minimum,maximum)
            elif (minimum != None):
                return self._check_min(x,minimum)
            else:
                return self._check_max(x,maximum)
        else:
            return x


    def multiply(self, factor, minimum = None, maximum = None, default_nodata_value = None):
        ''' 
        Multiply a grid with a factor or two grids with each other. Results in self grid.
        Arguments minimum and maximum are used to keep the result between minimum and maximum.
        '''
        if (type(factor) == Asciigrid):
            if (self.length == factor.length):
                # Check whether the is nodata in the grid
                if (self.__nodata_value == None and factor.__nodata_value == None):
                    self.values[:] = list(map(self._multi,self.values,factor.values))
                    # Check the range of the result
                    self._check_range_values(minimum,maximum)                    
                else:
                    # Check whether it is possible that nodata is allowed in self.grid
                    lwarn = False
                    lnodata_occur = False
                    if (self.nodata_value == None):
                        lwarn = True

                    for icell in range(self.length):
                        val1 = self.get_data(icell)
                        if (val1 != None):
                            val2 = factor.get_data(icell)
                            if (val2 != None):
                                val3 = self._check_range_value(val1*val2,minimum,maximum)
                                self.set_data(icell,val3)
                            else:
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the division by zero took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.MULTIPLY: Multiply could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)
                        else:
                            val2 = factor.get_data(icell)
                            if (val2 != None):
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the division by zero took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.MULTIPLY: Multiply could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)

                    # Check whether the result grid has nodata
                    if (lnodata_occur):
                        # Set grid properties for new nodata
                        self.nodata_value = default_nodata_value
                        self.nodata = True                      

            else:
                raise ASCIIGridError("ASCRASTER.MULTIPLY: Length of both grids is not the same.")

        elif (type(factor) == int or type(factor) == float):
            # Check whether the is nodata in the grid
            if (self.__nodata_value == None):
                self.values[:] = [factor*x for x in self.values]
                # Check the range of the result
                self._check_range_values(minimum,maximum)                    
            else:
                for icell in range(self.length):
                    val1 = self.get_data(icell)
                    if (val1 != None):
                        val2 = self._check_range_value(factor*val1,minimum,maximum)
                        self.set_data(icell,val2)
 
 
    def power(self, factor, minimum = None, maximum = None, default_nodata_value = None):
        ''' 
        Raise a grid to the power (factor) or two grids which calculates th epower. Results in self grid.
        Arguments minimum and maximum are used to keep the result between minimum and maximum.
        '''
        if (type(factor) == Asciigrid):
            if (self.length == factor.length):
                # Check whether the is nodata in the grid
                if (self.__nodata_value == None and factor.__nodata_value == None):
                    self.values[:] = list(map(self._power,self.values,factor.values))
                    if (len([x for x in self.values if x == None]) != 0):
                        # There are nodata values created in the calculation.
                        # Make nodata value for these grid cells.
                        if (default_nodata_value == None):
                            # Print also a traceback at this point to find to place where the power creates nodata took place
                            import traceback
                            print(traceback.print_exc())
                            raise ASCIIGridError("ASCRASTER.POWER: Power could result in None values. Specify default_nodata_value.")
                        else:
                            # Set grid properties for new nodata
                            self.nodata_value = default_nodata_value
                            self.nodata = True   
                    # Check the range of the result
                    self._check_range_values(minimum,maximum)                    
                else:
                    # Check whether it is possible that nodata is allowed in self.grid
                    lwarn = False
                    lnodata_occur = False
                    if (self.nodata_value == None):
                        lwarn = True
                
                    for icell in range(self.length):
                        val1 = self.get_data(icell)
                        if (val1 != None):
                            val2 = factor.get_data(icell)
                            if (val2 != None):
                                val3 = self._power(val1,val2)
                                if (val3 == None):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the power creates nodata took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.POWER: Power could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                val3 = self._check_range_value(val3,minimum,maximum)
                                self.set_data(icell,val3)
                            else:
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the power creates nodata took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.POWER: Power could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)
                        else:
                            val2 = factor.get_data(icell)
                            if (val2 != None):
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the power creates nodata took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.POWER: Power could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)
 
                    # Check whether the result grid has nodata
                    if (lnodata_occur):
                        # Set grid properties for new nodata
                        self.nodata_value = default_nodata_value
                        self.nodata = True                      
 
            else:
                raise ASCIIGridError("ASCRASTER.POWER: Length of both grids is not the same.")
 
        elif (type(factor) == int or type(factor) == float):
            # Check whether the is nodata in the grid
            if (self.__nodata_value == None):
                self.values[:] = [self._power(x,factor) for x in self.values]
                if (len([x for x in self.values if x == None]) != 0):
                    # There are nodata values created in the calculation.
                    # Make nodata value for these grid cells.
                    if (default_nodata_value == None):
                        # Print also a traceback at this point to find to place where the power creates nodata took place
                        import traceback
                        print(traceback.print_exc())
                        raise ASCIIGridError("ASCRASTER.POWER: Power could result in None values. Specify default_nodata_value.")
                    else:
                        # Set grid properties for new nodata
                        self.nodata_value = default_nodata_value
                        self.nodata = True 

                # Check the range of the result
                self._check_range_values(minimum,maximum)                    
            else:
                for icell in range(self.length):
                    val1 = self.get_data(icell)
                    if (val1 != None):
                        val2 = self._check_range_value(self._power(val1,factor),minimum,maximum)
                        self.set_data(icell,val2)
 
 
    def divide(self, factor, minimum = None, maximum = None, default_nodata_value = None):
        ''' 
        Divide a grid with a scalar or two grids with each other. Results (self / factor) in self grid.
        Arguments minimum and maximum are used to keep the result between minimum and maximum.
        '''
        if (type(factor) == Asciigrid):
            if (self.length == factor.length):
                # Check whether it is possible that nodata is allowed in self.grid
                lwarn = False
                lnodata_occur = False
                if (self.nodata_value == None):
                    lwarn = True

                for icell in range(self.length):
                    val1 = self.get_data(icell)
                    if (val1 != None):
                        val2 = factor.get_data(icell)
                        if (val2 != None):
                            try:
                                val3 = val1/val2
                                val3 = self._check_range_value(val3,minimum,maximum)
                                self.set_data(icell,val3)
                            except ZeroDivisionError:
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the division by zero took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.DIVIDE: Division could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)
                        else:
                            if (lwarn):
                                lnodata_occur = True
                                if (default_nodata_value == None):
                                    # Print also a traceback at this point to find to place where the division by zero took place
                                    import traceback
                                    print(traceback.print_exc())
                                    raise ASCIIGridError("ASCRASTER.DIVIDE: Division could result in None values. Specify default_nodata_value.")
                                else:
                                    val3 = default_nodata_value
                            else:
                                val3 = self.nodata_value
                            self.set_data(icell,val3)
                    else:
                        val2 = factor.get_data(icell)
                        if (val2 != None):
                            if (lwarn):
                                lnodata_occur = True
                                if (default_nodata_value == None):
                                    # Print also a traceback at this point to find to place where the division by zero took place
                                    import traceback
                                    print(traceback.print_exc())
                                    raise ASCIIGridError("ASCRASTER.DIVIDE: Division could result in None values. Specify default_nodata_value.")
                                else:
                                    val3 = default_nodata_value
                            else:
                                val3 = self.nodata_value
                            self.set_data(icell,val3)
                
                # Check whether the result grid has nodata
                if (lnodata_occur):
                    # Set grid properties for new nodata
                    self.nodata_value = default_nodata_value
                    self.nodata = True
            else:
                raise ASCIIGridError("ASCRASTER.DIVIDE: Length of both grids is not the same.")

        elif (type(factor) == int or type(factor) == float):
            try:
                dum = 1.0/factor
            except ZeroDivionError:
                raise ASCIIGridError("ASCRASTER.DIVIDE: Divide by zero is not allowed.")
            self.multiply(1.0/factor,minimum=minimum,maximum = maximum)


    def substract(self, factor, minimum = None, maximum = None, default_nodata_value = None):
        ''' 
        Substract a grid with a scalar or two grids with each other. Results (self - factor) in self grid.
        Arguments minimum and maximum are used to keep the result between minimum and maximum.
        '''
        if (type(factor) == Asciigrid):
            if (self.length == factor.length):
                # Check whether the is nodata in the grids
                if (self.__nodata_value == None and factor.__nodata_value == None):
                    self.values[:] = list(map(self._substraction,self.values,factor.values))
                    # Check the range of the result
                    self._check_range_values(minimum,maximum)                    
                else:
                    # Check whether it is possible that nodata is allowed in self.grid
                    lwarn = False
                    lnodata_occur = False
                    if (self.nodata_value == None):
                        lwarn = True

                    for icell in range(self.length):
                        val1 = self.get_data(icell)
                        if (val1 != None):
                            val2 = factor.get_data(icell)
                            if (val2 != None):
                                val3 = self._check_range_value(val1-val2,minimum,maximum)
                                self.set_data(icell,val3)
                            else:
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the division by zero took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.SUBSTRACTION: Substraction could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)
                        else:
                            val2 = factor.get_data(icell)
                            if (val2 != None):
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the division by zero took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.SUBSTRACTION: Substraction could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)

                    # Check whether the result grid has nodata
                    if (lnodata_occur):
                        # Set grid properties for new nodata
                        self.nodata_value = default_nodata_value
                        self.nodata = True

            else:
                raise ASCIIGridError("ASCRASTER.SUBSTRACT: Length of both grids is not the same.")

        elif (type(factor) == int or type(factor) == float):    
            self.add(-1.0 * factor,minimum=minimum,maximum = maximum)
            
    def maximum(self):
        '''
        Returns maximum value of this grid
        '''
        lfound = False
        maxvalue = None
        for icell in range(self.length):
            val = self.get_data(icell)
            if (val != None):
                if (lfound):
                    if (val > maxvalue):
                        maxvalue = val
                else:
                    maxvalue = val
                    lfound = True
        return maxvalue
      
    def minimum(self):
        '''
        Returns minimum value of this grid
        '''
        lfound = False
        minvalue = None
        for icell in range(self.length):
            val = self.get_data(icell)
            if (val != None):
                if (lfound):
                    if (val < minvalue):
                        minvalue = val
                else:
                    minvalue = val
                    lfound = True
        return minvalue
      
    def add(self, factor, minimum = None, maximum = None, default_nodata_value = None):
        ''' 
        Sum a grid with a scalar or two grids with each other. Results in self grid.
        Arguments minimum and maximum are used to keep the result between minimum and maximum.
        '''
        if (type(factor) == Asciigrid):
            if (self.length == factor.length):
                # Check whether the is nodata in the grids
                if (self.__nodata_value == None and factor.__nodata_value == None):
                    self.values[:] = list(map(self._summation,self.values,factor.values))
                    # Check the range of the result
                    self._check_range_values(minimum,maximum)                    
                else:
                    # Check whether it is possible that nodata is allowed in self.grid
                    lwarn = False
                    lnodata_occur = False
                    if (self.nodata_value == None):
                        lwarn = True

                    for icell in range(self.length):
                        val1 = self.get_data(icell)
                        if (val1 != None):
                            val2 = factor.get_data(icell)
                            if (val2 != None):
                                val3 = self._check_range_value(val1+val2,minimum,maximum)
                                self.set_data(icell,val3)
                            else:
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the division by zero took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.ADD: Addition could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)
                        else:
                            val2 = factor.get_data(icell)
                            if (val2 != None):
                                if (lwarn):
                                    lnodata_occur = True
                                    if (default_nodata_value == None):
                                        # Print also a traceback at this point to find to place where the division by zero took place
                                        import traceback
                                        print(traceback.print_exc())
                                        raise ASCIIGridError("ASCRASTER.ADD: Addition could result in None values. Specify default_nodata_value.")
                                    else:
                                        val3 = default_nodata_value
                                else:
                                    val3 = self.nodata_value
                                self.set_data(icell,val3)


                    # Check whether the result grid has nodata
                    if (lnodata_occur):
                        # Set grid properties for new nodata
                        self.nodata_value = default_nodata_value
                        self.nodata = True

            else:
                raise ASCIIGridError("ASCRASTER.ADD: Length of both grids is not the same.")

        elif (type(factor) == int or type(factor) == float):
            # Check whether the is nodata in the grid
            if (self.__nodata_value == None):
                self.values[:] = [factor+x for x in self.values]
                # Check the range of the result
                self._check_range_values(minimum,maximum)                    
            else:
                for icell in range(self.length):
                    val1 = self.get_data(icell)
                    if (val1 != None):
                        val2 = self._check_range_value(factor + val1,minimum,maximum)
                        self.set_data(icell,val2)


    def resize(self,number):
        '''
        Resize the grid to a bigger size.
        '''
        if (number < 1):
            raise ASCIIGridError("ASCRASTER.RESIZE: dimensions of new grid must be bigger.")
        if (self.mask != None):
            raise ASCIIGridError("ASCRASTER.RESIZE: No allowed for masked grid.")
        out = creategrid(self.ncols*number, self.nrows*number, self.xllcorner, self.yllcorner, \
                     self.cellsize/float(number), self.nodata_value,\
                     mask=None,numtype=self.numtype)
        for icell in range(out.length):
            irow,icol = out.get_row_col_from_index(icell)
            irow = irow//number
            icol = icol//number
            ind = self.get_index_from_row_col(irow,icol)
            out.set_data(icell,self.get_data(ind,self.nodata_value))
        return out
 
    def resample(self,number,method=1):
        '''
        Resample the smaller grid to a larger size. Sum the values for method = 1, method = 2 average.
        '''
        if (number < 1):
            raise ASCIIGridError("ASCRASTER.RESAMPLE: dimensions of new grid must be bigger.")
        if (self.mask != None):
            raise ASCIIGridError("ASCRASTER.RESAMPLE: No allowed for masked grid.")
        if (self.ncols%int(number) != 0):
            raise ASCIIGridError("ASCRASTER.RESAMPLE: Number must be a divider of the ncols.")
        if (self.nrows%int(number) != 0):
            raise ASCIIGridError("ASCRASTER.RESAMPLE: Number must be a divider of the nrows.")

        out = creategrid(self.ncols//number, self.nrows//number, self.xllcorner, self.yllcorner, \
                     self.cellsize*float(number), self.nodata_value,\
                     mask=None,numtype=self.numtype)
        if (method == 2):
            count = creategrid(self.ncols//number, self.nrows//number, self.xllcorner, self.yllcorner, \
                               self.cellsize*float(number), self.nodata_value,\
                               mask=None,numtype=self.numtype)

        for icell in range(self.length):
            val = self.get_data(icell)
            if (val != None):
                irow,icol = self.get_row_col_from_index(icell)
                irow = irow//number
                icol = icol//number
                ind = out.get_index_from_row_col(irow,icol)
                val_old = out.get_data(ind,0.0)
                out.set_data(ind,val_old+val)
                if (method == 2):
                    val_old = count.get_data(ind,0.0)
                    count.set_data(ind,val_old+1.)

        if (method == 2):
            for icell in range(out.length):
                try:
                    val = out.get_data(icell)
                    if (val != None):
                        out.set_data(ind,val/count.get_data(icell))
                except ZeroDivisionError:
                    pass

        return out


def compressed(filename):
    '''
    Checks whether the filename has an extension of .gz or .GZ. Then the file
    is compressed and the return value is True
    '''
    return os.path.splitext(filename)[1].upper() == ".GZ"

def read_header(filename):
    '''Read the header information from a ASCIIgrid.
    Returns a dictionary with the header properties.'''
    hps = {'ncols':None,
           'nrows':None,
           'xllcorner':None,
           'yllcorner':None,
           'cellsize':None,
           'nodata_value':None}

    # Opening grid ascii file
    try: 
        if (compressed(filename)):
            f = gzip.open(filename, 'r')
        else:      
            f = open(filename, 'r')
    except IOError:
        raise Exception("Error in opening file: %s" % filename)

    line = f.readline()
    # scan first 6 lines for parameters
    for _i in range(6):
        line_info = line.replace('\n', '').split()
        # Check if the found key is a valid keyword
        if line_info[0].lower() in list(hps.keys()):
            # set the value to the dictionary
            hps[line_info[0].lower()] = line_info[1]
        line = f.readline()
    # Check if all keywords and values are found and set
    for key in list(hps.keys()):
        # If one of the values is None, except nodata_value, the header was not completely read from file
        if hps[key] == None and key != 'nodata_value': raise ASCIIGridError("Incorrect header. Some attributes could not be found: %s" % key, f)
    f.close()
    return hps

def compare_grids(*args):
    '''
    Compares a given set of ascii raster files based on their headers.
    If the header properties are not equal False is returned.
    Checked properties are: ncols, nrows, xllcorner, yllcorner and cellsize.
    '''
    # Asume the given grids are equal.
    equal = True
    # If more than one raster is given compare them
    if len(args) > 1:
        # read the headers of the first one
        hps1 = read_header(args[0])
        # compare the others against the first
        hpsdic = {}
        # build a dictionary with properties
        for asciigrid in args[1:]:
            hpsdic[asciigrid] = read_header(asciigrid)
        for hpsn in list(hpsdic.values()):
            # if one header value deviates: set equal to False
            if hps1['ncols'] != hpsn['ncols']: equal = False
            if hps1['nrows'] != hpsn['nrows']: equal = False
            if hps1['xllcorner'] != hpsn['xllcorner']: equal = False
            if hps1['yllcorner'] != hpsn['yllcorner']: equal = False
            if hps1['cellsize'] != hpsn['cellsize']: equal = False
            # If a false 
            if not equal: break
        # print the found header info if the rasters are not equal
        if not equal:
            print("The given ascii rasters are not equal:")
            print("%s:" % args[0])
            print("ncols: %(ncols)s; nrows: %(nrows)s; xllcorner: %(xllcorner)s; yllcorner: %(yllcorner)s; cellsize: %(cellsize)s; " % hps1)
            for asciigrid in list(hpsdic.keys()):
                print("%s:" % asciigrid)
                print("ncols: %(ncols)s; nrows: %(nrows)s; xllcorner: %(xllcorner)s; yllcorner: %(yllcorner)s; cellsize: %(cellsize)s; " % hpsdic[asciigrid])
    return equal
    
def create_mask(ascii_file, value, logical = 'GE',numtype=None):
    '''
    Creates a mask. A mask is a list containing locators of all locations not masked.
    Before applying a mask on a Asciigrid instance compare the mask ascii_file to the instance.
    Returns a list.
    @param ascii_file: a valid filename to a ASCII raster file
    @type ascii_file: STRING
    @param value: Value of the cells not to be masked.
    @type value: FLOAT/INTEGER
    @param logical: Type of logical evaluation. Greater than or Equal is Default.
    GT - Greater than
    GE - Greater than or Equal
    LT - Less than
    LE - Less than or Equal
    EQ - Equal
    NQ - Not equal
    @type logical: STRING
    '''
    # TODO: make a mask class containing value_count to check if the mask can be used
    mask = []
    grid = Asciigrid(ascii_file=ascii_file,numtype=numtype)
    if logical == 'GT':
        for i in range(grid.length):
            if grid.values[i] > value:
                mask.append(i)
    elif logical == 'EQ':
        for i in range(grid.length):
            if grid.values[i] == value:
                mask.append(i)    
    elif logical == 'LT':
        for i in range(grid.length):
            if grid.values[i] < value:
                mask.append(i)
    elif logical == 'NQ':
        for i in range(grid.length):
            if grid.values[i] != value:
                mask.append(i)
    elif logical == 'GE':
        for i in range(grid.length):
            if grid.values[i] >= value:
                mask.append(i)
    elif logical == 'LE':
        for i in range(grid.length):
            if grid.values[i] <= value:
                mask.append(i)
    return mask

def readgrid(ascii_file, mask=None,numtype=None):
    '''
    Reads a ASCII raster from a file location. Returns a Asciigrid instance.
    @param ascii_file: a valid filename to a ASCII raster file
    @type ascii_file: STRING
    @param mask: a list containing pointers to non masked values 
    @type mask: LIST
    '''
    return Asciigrid(ascii_file=ascii_file, mask=mask,numtype=numtype)

def creategrid(ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value, mask=None,numtype=None):
    '''
    Creates a new ASCII raster based on the given parameters. Returns a Asciigrid instance.
    @param ncols: Number of columns
    @type ncols: INTEGER
    @param nrows: Number of rows
    @type nrows: INTEGER
    @param xllcorner: Lower left corner x coordinate
    @type xllcorner: FLOAT
    @param yllcorner: Lower left corner y coordinate
    @type yllcorner: FLOAT
    @param cellsize: Cellsize
    @type cellsize: FLOAT
    @param nodata_value: No data value
    @type nodata_value: INTEGER
    '''
    return Asciigrid(ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value, mask = mask,numtype=numtype)

def duplicategrid(asciigrid):
    '''
    Duplicates a Asciigrid instance.
    @param asciigrid: The Asciigrid instance to be copied
    @type asciigrid: Asciigrid
    '''
    dgrid = creategrid(asciigrid.ncols,
               asciigrid.nrows,
               asciigrid.xllcorner,
               asciigrid.yllcorner,
               asciigrid.cellsize,
               asciigrid.nodata_value,
               mask = asciigrid.mask,
               numtype=asciigrid.numtype)
    # Copy values into grid.
    dgrid.add_values(asciigrid.values)
      
    return dgrid
  
def interpolate_grid(time,time_start,grid_start,time_end,grid_end,nodata_interpol=0.0):
    '''
    Interpolation on time between two grids which are for time_start and time_end.
    When one of the grids has a nodata value than the nodata_interpol is used to perform a interpolation.
    This routine returns a ascgrid object with the interpolated values.
    If time is below time_start, then grid_start is returned, 
    if time is above time_end, then grid_end is returned.
    ''' 
    # Check given times.
    if (time_end < time_start):
        print("INTERPOLATE_GRID: given interpolation points is not correct.",\
              "First point must be smaller than second point.",\
              "First point found: " + str(time_start),\
              "Second point found: " + str(time_end))
        raise ASCIIGridError("INTERPOLATE_GRID: given interpolation points is not correct.")

    if (grid_start.length != grid_end.length):
        raise ASCIIGridError("INTERPOLATE_GRID: grids don't have the same length. Mask different?.")

    if (time <= time_start):
        return duplicategrid(grid_start)
    elif (time >= time_end):
        return duplicategrid(grid_end)
    else:
        # Interpolation
        factor = float(time - time_start)/float(time_end-time_start)
        grid_interpol = duplicategrid(grid_start) 
        for icell in range(grid_start.length):
            val_start = grid_start.get_data(icell)         
            val_end   = grid_end.get_data(icell)
            try:
                val_interpol = val_start + factor * (val_end - val_start)
            except TypeError:
                if (val_start == None and val_end == None):
                    # Do nothing
                    continue
                elif (val_start == None):
                    val_start = nodata_interpol
                else:
                    val_end = nodata_interpol
                val_interpol = val_start + factor * (val_end - val_start)
            # put interpolated value in output grid
            grid_interpol.set_data(icell,val_interpol)
        # Returns interpolated grid
        return grid_interpol

def checkCellSizesAreEqual(ListOfFileNames, log):
    
    #Check cellsize of all ascii-grids provided:
    print("Now starting CellSizeCheck")
    rasterCellSize = []
    for filename in ListOfFileNames:
        hps = read_header(filename)
        #print filename + "..."+ str(hps['cellsize'])
        rasterCellSize.append(hps['cellsize'])

    uniqueCellSizeValues = len(set(rasterCellSize))
    # print str(uniqueCellSizeValues)

    # In case all raster files have same cellsize: length of set is 1 (only one cellsize value present)
    if (uniqueCellSizeValues > 1):
        for filename in ListOfFileNames:
            hps = read_header(filename)
            print(filename + "..."+ str(hps['cellsize']))
            rasterCellSize.append(hps['cellsize'])

        log.write("ERROR ...Cellsizes do not match")
        print("ERROR ...Cellsizes do not match")
        #raise optparse.OptionValueError("Cellsizes do not match")
        return(0)
    "OK ...Cellsizes do match"
    return(1)

