# ******************************************************
## Revision "$LastChangedDate: 2018-07-13 17:07:05 +0200 (vr, 13 jul 2018) $"
## Date "$LastChangedRevision: 4 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/mocsy_module.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
'''
imports mocsy module
'''
#try:
#  import mocsy
#except ImportError:
import os
import sys
root = os.getenv("DGNM_ROOT")
if  root is None :
    print("***** ERROR *****")
    print("Environment parameter DGNM_ROOT is not set.")
    sys.exit(1)
if not os.path.isdir(root):
    print("***** ERROR ******")
    print("Environment parameter DGNM_ROOT is not set correctly.")
    print("Environment parameter DGNM_ROOT found: ",root)
    sys.exit(1)
# put mocsy directory in the sys.path.
if "/User" in root:
    print('Working in Mac OS')
    path = os.path.join(root,"libs","mocsy","mac")
    print('Fuck me')
elif "/theia" in root:
    print('Working in linux OS')
    path = os.path.join(root,"libs","mocsy","linux")
elif "/" in root:
    print('Working in linux OS - original statement')
    path = os.path.join(root,"libs","mocsy","linux")
if os.path.exists(path):
    sys.path.insert(3, path)
    