# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

#print("Import general_path.")
import os
import sys

# Read the name of the user for the personal version of modules.
name = 'carbon'
if (name == None):
    print("***** ERROR ******")
    print("Environment parameter DGNM_USER is not set.")
    sys.exit(1)
#print("NAME: ",name)

# Read the environment parameter with the rootdirectory
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
if (root == None):
    print("***** ERROR ******")
    print("Environment parameter DGNM_ROOT is not set.")
    sys.exit(1)
if (not os.path.isdir(root)):
    print("***** ERROR ******")
    print("Environment parameter DGNM_ROOT is not set correctly.")
    print("Environment parameter DGNM_ROOT found: ",root)
    sys.exit(1)
#print("ROOT: ",root)

# Read the environment parameter with the generalcode directory
generalcode = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'generalcode', 'trunk')
if (generalcode == None):
    print("***** ERROR ******")
    print("Environment parameter DGNM_GENERALCODE is not set.")
    sys.exit(1)
if (not os.path.isdir(generalcode)):
    print("***** ERROR ******")
    print("Environment parameter DGNM_GENERALCODE is not set correctly.")
    print("Environment parameter DGNM_GENERALCODE found: ",generalcode)
    sys.exit(1)
#print("GENERALCODE: ",generalcode)

# Set the generalcode directory in the python path
path = generalcode
if os.path.exists(path):
    sys.path.insert(0, path)
    print(path + " is added to the python search path for modules.") 

# Set the core directory in the python path
path = os.path.join(root,"core")
if os.path.exists(path):
    sys.path.insert(0, path)
    print(path + " is added to the python search path for modules.") 

# Set the personal directory in the python path
path = os.path.join(root,name,'code')
if os.path.exists(path):
    sys.path.insert(0, path)
    print(path + " is added to the python search path for modules.")
