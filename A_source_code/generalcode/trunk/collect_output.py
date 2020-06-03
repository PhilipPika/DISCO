# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/collect_output.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import sys
import os
import traceback

# Import Global ("generalcode") modules
# Generic code files ideally should be placed in the directory "generalcode"
# subfolder, but they may also be placed at the same base folder
__general = os.path.join(os.getcwd(), '/home/arthur/globalnutrients/generalcode/trunk')
#__general = os.path.join(r"../../..", 'generalcode/trunk')
if os.path.exists(__general): 
    sys.path.insert(0, __general)
    print(__general + " is added to the python search path for modules.") 
 
# Import general python modules
import optparse
import profile

# Import own general modules
#import aggregate
import general_class
import my_sys
import write_dict

outputfile = "overview"
scenarioname = "Anne"
aggregation_grid  = "continent9.map"
root_outputdirectory = os.getcwd()

total = {}

dirs = os.listdir(root_outputdirectory)
for idir in dirs:
    if (os.path.isdir(idir)):
        ldir = False
        # Check whether directory name is a year
        try:
            qq = int(idir)
            ldir = True
        except ValueError:
            pass

        if (ldir):
            print("Read year: ", idir)
            # Directory is a output directory which must be read.
            filename = os.path.join(idir,scenarioname +"_" + idir + "_" +  aggregation_grid + ".csv")
            if (os.path.isfile(filename)):
                # Read input of file
                total[idir] = general_class.read_general_file(filename,sep=";",out_type="list")
            else:
                 raise MyError("Filename: " + filename + " does not exist.")
            
# Sort all the years
keylist=[]
for key in total:
    keylist.append(key)
keylist.sort()

# Take the startyear of all the years
year_start = keylist[0]
# Take a name from the database
names_list = total[year_start][0].names

# Make a empty general class
G = general_class.General()
G.add_item("year",None)
for name in names_list:
    G.add_item(name,None)


# Make a list of all the output per aggregation region number
for iline in range(len(total[year_start])):
    # First time the header must be written.
    lheader = True
    # Write to output file
    fp = open(outputfile +  "_" + str(iline) + ".csv","w")
    for item in range(len(keylist)):
        out = G.copy()
        out = total[keylist[item]][iline].copy_into(out)
        out.set_val("year",keylist[item])
    
        out.write(fp,sep=";",lheader=lheader,NoneValue="")
        lheader = False
    
    fp.close()    
            




