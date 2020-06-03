# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/showimage.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import sys
import os
import ascraster
import my_sys

# Read input file
grid = ascraster.Asciigrid(ascii_file=sys.argv[1])
grid.save_as_bitmap("qq.jpeg",minV=0,maxV=5,mode="RGB",color=(1,0,0))
my_sys.my_system("display qq.jpeg")
my_sys.my_removefile("qq.jpeg")


