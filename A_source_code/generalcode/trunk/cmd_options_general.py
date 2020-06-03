# ******************************************************
## Revision "$LastChangedDate: 2019-02-01 08:45:22 +0100 (vr, 01 feb 2019) $"
## Date "$LastChangedRevision: 17 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/cmd_options_general.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import general python modules
import os
import sys
import optparse
import getpass
import time
import copy
import random

# Import own general modules
import ascraster
from error import *
import general_class

class Input(object):
    '''
    General functions for Parsing arguments and store them using the optparse module.

    Pass sys.argv into class: Input(sys.argv).
    Options will be included in 'options' object of instance of Input().
    Call options using:
    <instance>.options.<option>

    All other command-line arguments will be stored in <instance>.args
    '''

    # Switch debug messages on or off
    _debug = False

    def set_debug(self,val):
        self._debug = val

    def _dbg(self,msg):
        if (self._debug):
            print("DBG--> " + str(msg))

    def validate_file(self, file):
        '''
        Validate file

        Performs the following checks:
        * Check if file exists

        '''
        self._dbg("Checking file: %s" % file)
        if not os.path.exists(file):
            raise MyError("File '%s' does not exist" % file)

    def validate_directory(self, directory, bool_write=True):
        '''
        Validate directory.

        Performs the following checks:
        * Check if directory exists
        * Check if directory is really a directory
        * if bool_write is True Check for write access on directory
        and raise an exception if validation fails

        '''
        self._dbg("Checking directory: %s" % directory)
        if not os.path.exists(directory):
            raise MyError("Directory '%s' does not exist" % directory)
        elif not os.path.isdir(directory):
            raise MyError("Value '%s' is not a directory" % directory)
        elif bool_write:
            unique_jobid = "_" + str(int(random.random()*10000))
            test_filename = os.path.join(directory, "DeleteMe"+unique_jobid+".txt")
            try:
                file = open(test_filename,"w")
                file.close()
                os.remove(test_filename)
            except:
                raise MyError("You have no write access to '%s'" % directory)

    def _parse_inifile(self, list_args):
        '''
        Parse inifile and return list

        In inifile use command-line option notation, followed by '=', followed by argument, e.g. --inputdir = C:\temp.
        Comments must be preceded by '#'. Everything after '#' will be ignored by the parser.
        Command-line options prevail over inifile options.

        '''
        # Look for ini file in argument list
        inifile = list_args[list_args.index("--inifile")+1]
        # Check whether ini file exists
        self.validate_file(inifile)
        # Open ini file
        file = open(inifile,"r")
        ini_args_list = []
        for line in file.readlines():
            option = [elem.strip() for elem in line.split("#")[0].split("=")]
            # Check whether the first element is an option (first character is a dash)
            if option[0].startswith("-"):
                # If debug option is used in the ini file, activate the debug function
                if (option[0] == '-d' or option[0] == '--ldebug' and len(option) == 1):
                    self.set_debug(True)
                    option.append("1")
                elif (option[0] == '-d' or option[0] == '--ldebug'):
                    self.set_debug(self.make_boolean_value(option[1]))

                # Check the None type or logical value:
                elif (option[1] == "None"):
                    option[1] = None
                elif (option[1] == "True"):
                    option[1] = True
                elif (option[1] == "False"):
                    option[1] = False
                elif ("${" in option[1]):
                    # There is an environment variable used. Get the value.
                    option[1] = self.replace_environ_parameter(option[0],option[1])
                # Add to the list
                ini_args_list.extend(option)
        # Close ini file
        file.close()

        # If debug option is found, write the options found in the ini file to screen.
        for option in ini_args_list:
            self._dbg("Options or values found in inifile: %s" % option)

        return ini_args_list

    def _parse_parameter_inifile(self,cmd_class,inifile):
        '''
        Parse inifile with default parameters and return the class with this information

        In inifile use command-line option notation, followed by '=', followed by argument, e.g. P_bal_grass = "balance_P_grs.asc".
        Comments must be preceded by '#'. Everything after '#' will be ignored by the parser.
        The parameters are not set by the command-line options.

        '''
        # Check whether ini file exists
        self.validate_file(inifile)
        # Open ini file
        file = open(inifile,"r")
        for line in file.readlines():
            option = [elem.strip() for elem in line.split("#")[0].split("=")]
            if (len(option) == 2):
                # Add to the cmd_class
                try:
                    val = int(option[1])
                except ValueError:
                    try:
                        val = float(option[1])
                    except ValueError:
                        val = option[1]
                        # Check the None type or logical value:
                        if (option[1] == "None"):
                            val = None
                        elif (option[1] == "True"):
                            val = True
                        elif (option[1] == "False"):
                            val = False
                        elif ("${" in option[1]):
                            # There is an environment variable used. Get the value.
                            val = self.replace_environ_parameter(option[0],option[1])
                setattr(cmd_class,option[0],val)
        
        # Close ini file
        file.close()

    def replace_environ_parameter(self,varname,val):
        '''
        There is an environment variable used. Get the value and return the value with the value 
        of the environment parameter.
        If environment parameter is not found, an error will occurr.
        '''
        ibegin = val.index("{")
        iend = val.index("}")
        envval = os.getenv(val[ibegin+1:iend])
        #print("CMD line: ",varname,envval,val[ibegin:iend+1],val)
        if (envval != None):
            outname = val[:ibegin-1] + envval + val[iend+1:]
            #print("CMD line: ",varname,outname)
            lfound = (outname.find("{") != -1)
            if (lfound):
                outname = self.replace_environ_parameter(varname,outname)
            return outname
        else:
            raise MyError("Environment parameter " + str(val[ibegin+1:iend]) + " is not found.",\
                          "This is used in model parameter: " + varname)
        lfound = (val.find("{") != -1)
        if (lfound):
           val = self.replace_environ_parameter(varname,val)

    # Conversion to a boolean for local stuff:
    def make_boolean_value(self,val):
        '''
        Make from value a boolean
        '''
        try:
            val = int(val)
            return val
        except ValueError:
            # Check the None type or logical value:
            if (val == "None"):
                return False
            elif (val == "True"):
                 return True
            elif (val == "False"):
                 return False
            else:
                 return False

    # Conversion to a boolean:
    def make_boolean(self,name,lnot_present=False):
        try:
            a = getattr(self.options,str(name))
            if (type(a) == str):
                raise MyError(name + " must be a boolean.")
            elif (type(a) != bool):
                if (type(a) != int):
                    setattr(self.options,str(name),False)
                elif (self.options.lcompress != 1):
                    setattr(self.options,str(name),False)
                else:
                    setattr(self.options,str(name),True)
        except AttributeError:
            setattr(self.options,str(name),lnot_present)
 
    
    # Following functions are used to perform a sensitivity analyse.
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
        # Check the range of the result
        if (val_min != None or val_max != None):
            if (val_min != None and val_max != None):
                x1 = self._check_min(x,val_min)
                return self._check_max(x1,val_max) 
            elif (val_min != None):
                return self._check_min(x,val_min)
            else:
                return self._check_max(x,val_max)
        else:
             return x 
                     
    def label_present(self,name):
        try:          
            a = getattr(self.options,str(name))
            return True
        except AttributeError:
            return False

    def grid_mult_scalar(self,label,filename,xmin=None,xmax=None):
        '''
        Make a new grid file for a sensitivity analyse. Filename is changed to output directory.
        '''
        if (self.label_present(label)):
            print(label + ' is found.')
            # Store old filename
            filename_old  = filename            
            filename  = os.path.join(self.options.outputdir,os.path.basename(filename_old))
            # Read input file
            grid=ascraster.Asciigrid(ascii_file=filename_old)
            grid.multiply(float(getattr(self.options,label)),minimum=xmin,maximum=xmax)
            grid.write_ascii_file(filename)
            del grid
        else:
            print(label + ' is NOT found.')

        # Return new file name
        return filename        

    def grid_add_scalar(self,label,filename,minimum=None, maximum=None):
        '''
        Make a new grid file for a sensitivity analyse. Filename is changed to output directory.
        '''
        if (self.label_present(label)):
            print(label + ' is found.')
            # Store old filename
            filename_old  = filename
            filename  = os.path.join(self.options.outputdir,os.path.basename(filename_old))
            # Read input file
            grid=ascraster.Asciigrid(ascii_file=filename_old)
            grid.add(float(getattr(self.options,label)),minimum=minimum,maximum=maximum)
            grid.write_ascii_file(filename)
            del grid
        else:
            print(label + ' is NOT found.')

        # Return new file name
        return filename

    def table_mult_scalar(self,label,filename,key="isocode",sep=";",minimum=None, maximum=None):
        '''
        Make a new table file for a sensitivity analyse. Filename is changed to output directory.
        '''
        if (self.label_present(label)):
            print(label + ' is found.')
            # Store old filename
            filename_old  = filename
            filename  = os.path.join(self.options.outputdir,os.path.basename(filename_old))
            # Read input file
            data_dict = general_class.read_general_file(filename_old,\
                                                sep=sep,out_type="list")

            # Get multiplier
            mult_factor = float(getattr(self.options,label))
            
            # Multiply every element of the file
            for item in range(len(data_dict)):
                for name in data_dict[item].names:
                    if (name == key):
                        continue
                    try:
                        val = float(data_dict[item].get_val(name))
                        #print data_dict[item].get_val("isocode"),name,val,
                        val = self._check_range(mult_factor*val,minimum,maximum)
                        #print val
                        data_dict[item].set_val(name,val)
                    except ValueError as TypeError:
                        pass

            # Write all info to file
            fp = open(filename,"w")
            # Write header to file
            fp.write(sep.join(data_dict[0].names)+"\n")

            # Write data block to file.
            lastname = data_dict[0].names[-1]
            for item in range(len(data_dict)):
                for name in data_dict[item].names[:-1]:
                    fp.write(str(data_dict[item].get_val(name))+sep)
                fp.write(str(data_dict[item].get_val(lastname))+"\n")
            fp.close()

        else:
            print(label + ' is NOT found.')        

        # Return new file name
        return filename

def make_parameter_list(obj,basename,lfatal=True):
    '''
    This function puts the given information from the parameters input file into one structure.
    All parameters of params.<basename>* are joined together in a list.
    The joined information is put back in the params.<basename>_list object.
    '''
    try:
        # Check whether this parameter is already set. 
        # When it is set, then do nothing.
        #print str(basename)+"_list"
        qq = getattr(obj,str(basename)+"_list")
    except AttributeError:
        counter = 1
        listout = []
        lfound = True
        while (lfound):
            try:
                #print basename + str(counter),getattr(obj,str(basename)+str(counter))
                val = getattr(obj,str(basename)+str(counter))
                listout.append(val)
                counter += 1
            except AttributeError:
                lfound = False

        # Check whether there is information found.
        if (listout == []):
            if (lfatal):
                raise MyError("There is no information found like "+str(basename)+"1 in the parameters object.")

        # Put new list in the params object    
        setattr(obj,str(basename)+"_list",copy.deepcopy(listout))
        

def make_weighing_list(obj,basename):
    '''
    Get all information of the weighing maps and properties and put this in the class obj.
    '''
    # Read the maps
    make_parameter_list(obj,basename,lfatal=False)
    
    # Get for each map the time values and the weighing values.
    for item in range(len(getattr(obj,basename+"_list"))):
        make_parameter_list(obj,basename+str(item+1)+"_t")
        make_parameter_list(obj,basename+str(item+1)+"_y")
        
        # Check input
        if (len(getattr(obj,basename+str(item+1)+"_t_list"))  != len(getattr(obj,basename+str(item+1)+"_y_list"))):
            raise MyError("Information to compile "+basename+" in the parameters file is not correct.",\
                          "Number of entries for "+basename+str(item+1)+"_t_list = "+str(len(getattr(obj,basename+str(item+1)+"_t_list"))),\
                          "Number of entries for "+basename+str(item+1)+"_y_list = "+str(len(getattr(obj,basename+str(item+1)+"_y_list")))) 

        
def make_ref_list(obj,basename):
    '''
    Get all information of the reference maps and properties and put this in the class obj.
    '''
    # Read the maps
    make_parameter_list(obj,basename,lfatal=True)

    # Get for each map the time values and the weighing values.
    make_parameter_list(obj,basename+"_t")
    
    # Check input
    if (len(getattr(obj,basename+"_list"))  != len(getattr(obj,basename+"_t_list"))):
        raise MyError("Information to compile "+basename+" in the parameters file is not correct.",\
                          "Number of entries for "+basename+"_t_list = "+str(len(getattr(obj,basename+"_t_list"))),\
                          "Number of entries for "+basename+"_list = "+str(len(getattr(obj,basename+"_list")))) 
