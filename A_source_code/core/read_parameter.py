# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/read_parameter.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import Python modules
import os

# Import own general modules
import general_class
from error import *

# Import local modules.
import specie
import general_func

def readfile(filename):
    '''
    Read all the specs for the species and processes.
    Returns a list of species and a list of processes.
    '''
    proclist = []
    srclist = [] #added LV 22-07-2016: description of load sources in inifile 
    speclist = []
    paramlist = []
    if not os.path.exists(filename):
        raise MyError("File '%s' does not exist" % filename)

    # Open file with all parameters.
    fp = open(filename,"r")

    # Determine all the information per line.
    obj = None
    for line in fp.readlines():
        option = [elem.strip() for elem in line.split("#")[0].split("=")]
        if (len(option[0]) < 1):
            # Empty line
            continue

        # Check whether the first element is start of a new class 
        # (first character is a "[" ). 
        if option[0].startswith("["):
           # Start of a new class
           if (option[0].strip().upper() == "[PROCES]"):
               # Save old obj
               if (obj != None):
                   if (obj.get_val("type") == "GENERAL"):
                       proclist.append(convert(obj))
                   elif (obj.get_val("type") == "SRC"):
                       srclist.append(convert(obj))
                   else:
                       speclist.append(convert(obj))
               obj = general_class.General()
               obj.add_item("type","GENERAL")
           elif (option[0].strip().upper() == "[SOURCE]"): #added LV 22-07-2016
               # Save old obj
               if (obj != None):
                   if (obj.get_val("type") == "GENERAL"):
                       proclist.append(convert(obj))
                   elif (obj.get_val("type") == "SRC"):
                       srclist.append(convert(obj))
                   else:
                       speclist.append(convert(obj))
               obj = general_class.General()
               obj.add_item("type","SRC")
           elif (option[0].strip().upper() == "[SPECIE]"):
               # Save old obj
               if (obj != None):
                   if (obj.get_val("type") == "GENERAL"):                   
                       proclist.append(convert(obj))
                   elif (obj.get_val("type") == "SRC"):
                       srclist.append(convert(obj))
                   else:                   
                       speclist.append(convert(obj))                       
               obj = general_class.General()
               obj.add_item("type","SPEC")                                                          
           elif (option[0].strip().upper() == "[PARAMETERS]"):
               # Save old obj
               if (obj != None):
                   if (obj.get_val("type") == "GENERAL"):
                       proclist.append(convert(obj))
                   elif (obj.get_val("type") == "SRC"):
                       srclist.append(convert(obj))
                   else:
                       speclist.append(convert(obj))
               obj = general_class.General()
               obj.add_item("type","GENERAL")
               # Set extra item, so this entries can easily extracted out of proclist
               obj.add_item("subtype","PARAMETERS")
                              
           else:
               raise MyError("File '%s' contains errors. Class %s is not known." % (filename,option[0]))

        else:
            # Here information from a class is given.
            try:
                obj.add_item(option[0],option[1])
            except IndexError:
                raise MyError("File '%s' contains errors." % filename,\
                              "Can not handle line: %s" % line)
    # Close input file
    fp.close()

    # Add last class to output list
    if (obj != None):
        if (obj.get_val("type") == "GENERAL"):
            proclist.append(convert(obj))
        elif (obj.get_val("type") == "SRC"):
            srclist.append(convert(obj))
        else:
            speclist.append(convert(obj))            
            
    # Check whether the processes are refering to the right species and extract the general parameters.
    remove_list = []
    for iproc in range(len(proclist)):
        obj = proclist[iproc]
        try:
            qq = obj.get_val("subtype")
            # This is a parameter entry
            remove_list.append(iproc)
            continue
        except:
            pass
        for name in obj.get_attrib():
            if (name.startswith("dy_")):
                try:
                    ind = general_func.find_spec(speclist,name[3:])
                except MyError:
                    raise MyError("Process %s refers to specie which does not exist." % (obj.get_val("name")),\
                                  "Speciesname found: %s" % name[3:])

    # Check whether the sources are refering to the right species
    for isrc in range(len(srclist)):
        obj = srclist[isrc]
        for name in obj.get_attrib():
            if (name.startswith("fr_")):
                try:
                    ind = general_func.find_spec(speclist,name[3:])
                except MyError:
                    raise MyError("Source %s refers to specie that does not exist." % (obj.get_val("name")),\
                                  "Speciesname found: %s" % name[3:])
                
    # Make one general object for the parameters.
    params_local = general_class.General()
    for item in remove_list:
        params_local.merge(proclist[item])
    # Correct the names list of general class
    params_local.names = list(set(params_local.names))
    
    # Delete the parameters entries in the proclist
    for item in range(len(remove_list)-1,-1,-1):
        del proclist[remove_list[item]]
        
    #return speclist,proclist,params_local
    return speclist,srclist,proclist,params_local
    
def convert(obj):
    '''
    Conversion from General() to another type of class.
    '''
    try:
        if (obj.get_val("type") == "GENERAL"):
            # It is a General class.
            outobj = general_class.General()
        elif (obj.get_val("type") == "SRC"):
            # It is a General class.
            outobj = general_class.General()
        elif (obj.get_val("type") == "SPEC"):
            # It is a Spec class.
            outobj = specie.Spec()                  
        else:
            raise MyError("It is not possible to convert class of type '%s' . Class is not known." % obj.get_val("type"))
    except AttributeError:
        raise MyError("It is not possible to convert this object.")

    # Put all attributes in the output class
    for name in obj.get_attrib():
        if (obj.get_val("type") == "GENERAL") or (obj.get_val("type") == "SRC"):
            outobj.add_item(name,str2num(obj.get_val(name)))
        else:
            try:
                # Check whether this attribute is part of this class.
                dum = getattr(outobj,name)                
            except AttributeError:
                raise MyError("Attribute %s is unknown for class of type: %s." % (name,obj.get_val("type")))
            setattr(outobj,name,str2num(obj.get_val(name)))      
    try:
        # Change all attributes from strings into the type that is required.
        #outobj.check_type()
        pass
    except AttributeError:
        raise MyError("Method check_type is not found for class of type: %s." % obj.get_val("type"))

    return outobj

def str2num(s):
    '''
    checks if a string is a float or integer. Returns a float or integer or a string.
    '''
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s
