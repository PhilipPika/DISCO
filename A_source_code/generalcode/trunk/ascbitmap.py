# ******************************************************
## Revision "$LastChangedDate: 2019-01-31 15:27:30 +0100 (do, 31 jan 2019) $"
## Date "$LastChangedRevision: 16 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/ascbitmap.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
'''
ascbitmap.py

A Python library containing functions to read and write bitmap files.

Initialy created on 26 mei 2009
@author: warrinka
'''
try:
    from PIL import Image
    PILavailable = True
except ImportError:
    PILavailable = False

# General python modules
import os
import math

# Generalcode modules
from iround import *
from error import *

if PILavailable:
    class Asciibitmap:
        '''
        Creates a graphic representation of an asciigrid.
        '''
        def __init__(self,mode=None):
            self.nrows = None
            self.ncols = None
            self.values = None
            self.nodata = None
            # Set which colorband (gray, RGB or CMYK) is used.
            if (mode == None):
                self.mode = "L"
            else:
                self.mode = mode.upper()
            self.nodatacolor = None

        def createbitmap(self, nrows, ncols, values, filename, nodata = None, minV = None, maxV = None,\
                         color=None, number_of_classes = None, classtype = None, nodatacolor=None):
            '''
            Creates a bitmap-file with a graphic representation of the given values.
            @param nrows: Number of rows
            @type nrows: INTEGER
            @param ncols: Number of columns
            @type ncols: INTEGER
            @param values: List with values (INTEGER or FLOAT)
            @type values: LIST
            @param filename: Path including filename for output file
            @type filename: STRING
            @param nodata: (Optional) No_data value.
            @type nodata: FLOAT/INTEGER
            @param color: (Optional) Tuple with active colors.
            @type nodata: Tuple of three integers. (red, green, blue)
            @param number_of_classes: (Optional) Number of classes
            @type number_of_classes: INTEGER
            @param classtype: (Optional) Classtype specification
            @type classtype: Text with LOG for logaritmic distribution of the classes.
                             Text with filename, with on each line value and colorcode.
                             The value is max boundary (boundary included) for that color. 
                             (in case of RGB three integers otherwise one integer) 
            @nodatacolor: Color for the nodata
            @type nodatacolor: Tuple of three integers (in RBG mode) 
            '''
            lfixed_class = False
            self.nrows = nrows
            self.ncols = ncols
            self.values = values
            self.nodata = nodata
            self.nodatacolor = nodatacolor
            self.size = (self.ncols, self.nrows)
            self.img = Image.new(self.mode, self.size)
            if (maxV == None):
                self.maxV = max(list(filter(self.isno_nodata,self.values)))
            else:
                self.maxV = maxV
            if (minV == None):
                self.minV = min(list(filter(self.isno_nodata,self.values)))
            else:
                self.minV = minV
            
            if (classtype != None):
                if (classtype.upper() == "LOG"):
                    # Shift everything by one to avoid zero's.
                    self.minV = math.log(self.minV+1)
                    self.maxV = math.log(self.maxV+1)
                    classtype = classtype.upper()
                else:
                    # Classtype is a filename
                    # Number of classes is not used in this situation.
                    if (os.path.isfile(classtype)):
                        lfixed_class = True 
                        fp = open(classtype,"r")
                        lines = fp.readlines()
                        fp.close()
                        boundary = []
                        color_class = []
                        for item in range(len(lines)):
                            fields = lines[item].split()
                            if (len(fields) == 0):
                                continue
                            try:
                                if (self.mode == "L"):
                                    boundary.append(float(fields[0]))
                                    color_class.append(int(fields[1]))
                                elif (self.mode == "RGB"):
                                    boundary.append(float(fields[0]))
                                    color_class.append((int(fields[1]),int(fields[2]),int(fields[3])))
                            except IndexError:
                                print("***** ERROR *****")
                                print("File " + classtype + " is not correct.")
                                print("Last line read: ",lines[item]) 
                                raise IOError
                            
                    else:
                        print("***** ERROR *****")
                        print("File " + classtype + " is not found.") 
                        raise IOError
                
            if (lfixed_class):
                for r in range(self.nrows):
                    for c in range(self.ncols):
                        val_index = (r*self.ncols)+c
                        if self.values[val_index] == nodata or self.values[val_index] == None:
                            if (self.nodatacolor != None):
                                value = self.nodatacolor
                            else:
                                if (self.mode == "L"):
                                    value = 0
                                elif (self.mode == "RGB"):
                                    value = (0,0,0)
                        else:
                            for iclass in  range(len(boundary)):
                                if (self.values[val_index] <= boundary[iclass]):
                                    value = color_class[iclass]
                                    break
                            else:
                                # Value is bigger than the last given number.
                                if (self.mode == "L"):
                                    value = 255
                                elif (self.mode == "RGB"):
                                    #value = (255,255,255)
                                    value = color_class[-1] #LV add an option somewhere...

                        self.img.putpixel((c,r), value)           
               
            elif (self.mode == "L"):
                for r in range(self.nrows):
                    for c in range(self.ncols):
                        val_index = (r*self.ncols)+c
                        if self.values[val_index] == nodata or self.values[val_index] == None:
                            if (self.nodatacolor != None):
                                value = self.nodatacolor
                            else:
                                value = 0
                        else:
                            if (classtype == "LOG"):
                                val = math.log(self.values[val_index]+1)

                            else:
                                val = self.values[val_index]
                            value = self._get_value_255(self.maxV, self.minV, val,\
                                                        number_of_classes = number_of_classes)

                        self.img.putpixel((c,r), value)
            elif (self.mode == "RGB"):
                if (int(color[0]) == 1):
                    red = 1
                else:
                    red = None
                if (int(color[1]) == 1):
                    green = 1
                else:
                    green = None
                if (int(color[2]) == 1):
                    blue = 1
                else:
                    blue = None

                for r in range(self.nrows):
                    for c in range(self.ncols):
                        val_index = (r*self.ncols)+c
                        if self.values[val_index] == nodata or self.values[val_index] == None:
                            if (self.nodatacolor != None):
                                value = self.nodatacolor
                            else:
                                value = (0,0,0)
                        else: 
                            if (classtype == "LOG"):
                                val = math.log(self.values[val_index]+1)
                            else:
                                val = self.values[val_index]
                            value = self._get_value_RGB(self.maxV, self.minV, val,\
                                                        red=red,blue=blue,green=green,\
                                                        number_of_classes = number_of_classes)
                        self.img.putpixel((c,r), value)
            
            self._savebitmap(filename)


        def _savebitmap(self, filename):
            '''
            Saves the image under the given filename. The format is determined from the filename extension.
            @param filename: Path including filename for output file
            @type filename: STRING
            '''
            self.img.save(filename)

        def _get_value_RGB(self, maxV, minV, val,red=None,blue=None,green=None,\
                           number_of_classes = None):
            '''
            Scales a given value on a scale from 0 - 255
            @param maxV: Maximum value from values in data
            @type maxV: FLOAT/INTEGER
            @param minV: Minimum value from values in data
            @type minV: FLOAT/INTEGER
            @param val: Value 
            @type val: FLOAT/INTEGER
            @param red: a value. With a value, the red range is used.
            @type red: FLOAT/INTEGER
            @param blue: a value. With a value, the blue range is used. 
            @type blue: FLOAT/INTEGER
            @param green: a value. With a value, the green range is used. 
            @type green: FLOAT/INTEGER
            '''
            dif = maxV - minV
            if (number_of_classes == None):
                try:
                    step = 200.0/dif
                except ZeroDivisionError:
                    step = 1
            else:
                try:
                    step = float(number_of_classes)/dif
                except ZeroDivisionError:
                    step = 1

            if (red != None):
                if (blue != None):
                    if (green != None):
                        # Go from blue to green to red
                        if (val - minV > 0.5*dif):
                            value = (min(max(0,int((val-0.5*dif-minV)*0.5*step)),200)+55,\
                                    min(max(0,int(200 -(val-0.5*dif-minV)*0.5*step)),200)+55,\
                                    0)
                        else:
                            value = (0,\
                                    min(max(0,int((val-0.5*dif-minV)*0.5*step)),200)+55,\
                                    min(max(0,int(200 -(val-0.5*dif-minV)*0.5*step)),200)+55)
                    else:
                        # For blue to red
                        value = (self.interpol_up(step,minV,val),\
                                 0,\
                                 self.interpol_down(step,minV,val))
              
                else:
                    # It is red, not blue
                    if (green != None):
                       # For blue to red
                        value = (self.interpol_up(step,minV,val),\
                                 self.interpol_down(step,minV,val),\
                                 0)
                    else:
                       # Only red (from black to red).
                        value = (self.interpol_up(step,minV,val),\
                                 0,\
                                 0)
            else:
                # No red
                if (blue != None):
                    if (green != None):
                        # Go from blue to green
                        value = (0,\
                                 self.interpol_up(step,minV,val),\
                                 self.interpol_down(step,minV,val))
                    else:
                        # For blue (from black to blue)
                        value = (0,\
                                 0,\
                                 self.interpol_up(step,minV,val))
              
                else:
                    # It is not red, not blue
                    if (green != None):
                       # For green (from black to green)
                        value = (0,\
                                 self.interpol_up(step,minV,val),\
                                 0)
                    else:
                       # Gray from light to dark.
                        value = (self.interpol_up(step,minV,val),\
                                 self.interpol_up(step,minV,val),\
                                 self.interpol_up(step,minV,val))


            return value

        def interpol_up(self,step,minV,val):
            return min(max(0,int((val-minV)*step)),200) + 55

        def interpol_down(self,step,minV,val):
            return min(max(0,int(200-(val-minV)*step)),200) + 55
            
        def _get_value_255(self, maxV, minV, val, number_of_classes = None):
            '''
            Scales a given value on a scale from 0 - 255
            @param maxV: Maximum value from values in data
            @type maxV: FLOAT/INTEGER
            @param minV: Minimum value from values in data
            @type minV: FLOAT/INTEGER
            @param val: Value
            @type val: FLOAT/INTEGER
            '''
            dif = maxV - minV
            if (number_of_classes == None):
                try:
                    step = 200.0/dif
                except ZeroDivisionError:
                    step = 1
            else:
                try:
                    step = float(number_of_classes)/dif
                except ZeroDivisionError:
                    step = 1

            value = (val - minV)*step
            return min(max(0,int(value)),200) + 55
        
        def isno_nodata(self,val):
            return (val != self.nodata)
        
        
        def showbitmap(self):
            '''
            Saves the image to a temporary BMP file, and uses the standard BMP display utility to show it.
            '''
            self.img.show()
else: 
    print("WARNING:\nPython Imaging Library (PIL) is not installed.\nAsciibitmap class will not be available.")


def make_legend_text(colorfile,starttext="nodata"):
    '''
    Read file with colorsettings and extract the legend text.
    '''
    legend_text = [starttext]
    
    # Open color settings file
    fp = open(colorfile,"r")
    lines = fp.readlines()
    fp.close()

    boundaries = []
    for line in lines:
        fields = line.split()
        if (len(fields) == 0):
            # Empty line, goto next
            continue
        try:
            boundaries.append(fields[0])
        except:
            import traceback
            traceback.print_exc()
            raise MyError("Reading legend text from file: " + str(colorfile) + " goes wrong.")

    for item in range(len(boundaries) - 2):
        legend_text.append(str(boundaries[item]) +" - " + str(boundaries[item+1]))

    # The last element is different
    legend_text.append("> " + str(boundaries[-2]))

    return legend_text

def make_legend_color(colorfile):
    '''
    Read file with colorsettings and extract the legend text.
    '''
    legend_color = []
    
    # Open color settings file
    fp = open(colorfile,"r")
    lines = fp.readlines()
    fp.close()

    for line in lines:
        fields = line.split()
        if (len(fields) == 0):
            # Empty line, goto next
            continue
        try:
            legend_color.append((iround(float(fields[1])),iround(float(fields[2])),iround(float(fields[3]))))
        except:
            import traceback
            traceback.print_exc()
            raise MyError("Reading legend colors from file: " + str(colorfile) + " goes wrong.")


    return legend_color

def row_conversion(val,cellsize_invers,yll,nrows):
    return nrows -(cellsize_invers*(val-yll))

def column_conversion(val,cellsize_invers,xll):
    return cellsize_invers*(val-xll)

def draw_polygons(im,grid,polygonfile,linecolor=(0,0,0)):
    '''
    Read a XY polygon file and draw this over the image.
    @im : Pointer to the image which will be changed.
    @grid: ascraster object of the raster that is the bitmap of (we only use the header of this grid).
    @polygonfile: File with on each line the label of the polygon or coordinates of one point of the polygon (XY format). 
    @linecolor: Tuple with the RGB colors.
    '''
    if (polygonfile == None):
        # Nothing to do.
        return
    elif not (os.path.isfile(polygonfile)):
        print("Warning: Polygonfile: " + polygonfile + " is not found.")
        return

    # Read the polygon file.
    fp = open(polygonfile,"r")
    lines = fp.readlines()
    fp.close()

    # Make drawing possible for this figure
    draw = ImageDraw.Draw(im)

    # Get info out of the headr of the grid.
    xll = float(grid.xllcorner)
    yll = float(grid.yllcorner)
    cellsize_invers = 1.0/float(grid.cellsize)
    nrows = iround(grid.nrows)

    # Draw each polygon
    polyset = []
    for line in lines:
        fields = line.split()
        if (len(fields)<=1):
            # We can draw this polygon.
            if (len(polyset)> 0):
                draw.polygon(polyset, outline=linecolor)            
            # Make polyset empty again
            polyset = []
        else:
            # Conversion of the coordinates and put in the list.
            polyset.append(tuple([column_conversion(float(fields[0]),cellsize_invers,xll),\
                                  row_conversion(float(fields[1]),cellsize_invers,yll,nrows)]))

    # Draw last polygon
    if (len(polyset)> 0):
        draw.polygon(polyset, outline=linecolor)            
        
    del draw 

def legend(im, upleft_pos, boxsize, Font,legend_color, legend_text,legend_title="",fontsize=40):
    '''Adds a time slider with current position to a image.
    
    @param im: Image instance to draw the legend on
    @type im: Image instance
    @param upleft_pos: Coordinate left upper corner of legend box
    @type upleft_pos: tuple/list
    @param width: Width of the legend box
    @type width: integer
    @param length: Length of the legend box
    @type length: integer
    @param legend_color: List with color numbers
    @type legend_color: list integers
    @param legend_text: List with legend text
    @type legend_text: list of strings
    '''
    draw = ImageDraw.Draw(im)
    # Size of the colorsquares:
    pos = upleft_pos
    font = ImageFont.truetype(Font, fontsize)
    if (legend_title != ""):
        # Write title to image
        draw.text((pos[0],pos[1]+boxsize-10),legend_title,(0,0,0),font=font)
        pos = (pos[0],pos[1]+boxsize)

    for item in range(len(legend_color)):
        pos = (pos[0],pos[1]+boxsize)
        # Draw colorbox
        draw.polygon([pos, (pos[0]+boxsize,pos[1]), (pos[0]+boxsize,pos[1]+boxsize),\
                     (pos[0],pos[1]+boxsize)],fill=legend_color[item], outline=(0,0,0))
        draw.text((pos[0]+boxsize+15,pos[1]+int(0.25*boxsize)),legend_text[item],(0,0,0),font=font)
    del draw
