# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/convert_grid_dir.py $"
# ******************************************************
'''
Test script to test the functionality of ascraster.
'''
import os
import sys
__general = r"/home/beusena/globalnutrients/generalcode/trunk"
if os.path.exists(__general): 
    sys.path.insert(0, __general)
    print(__general + " is added to the python search path for modules.") 

import ascraster
import my_sys

inputdir = sys.argv[1]
outputdir = sys.argv[2]


if not os.path.isdir(outputdir):
    os.makedirs(outputdir)

# Convert total directory of raster files to full ascii grid files.
listdir = os.listdir(inputdir)

for file in listdir:
    filename = os.path.join(inputdir,file)
    if (os.path.splitext(file)[1].upper() == ".ASC"):
        # Convert this file
        print(filename + " is converted.")       
        grid = ascraster.Asciigrid(ascii_file=filename)
        fileout = os.path.join(outputdir,file)
        grid.write_ascii_file(fileout)
    elif (os.path.splitext(file)[1].upper() == ".GZ"):
        if (os.path.splitext(os.path.splitext(file)[0])[1].upper() == ".ASC"):
            # Convert this file 
            print(filename + " is converted.")   
            grid = ascraster.Asciigrid(ascii_file=filename)
            fileout = os.path.join(outputdir,os.path.splitext(file)[0])
            grid.write_ascii_file(fileout)
        else:
            print(filename + " is skipped.") 
    else:
        # Copy this file to new directory
        #my_sys.my_copyfile(filename,fileout)
        print(filename + " is skipped.")
        pass


