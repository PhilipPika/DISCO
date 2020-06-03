# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# General python modules
import os
import sys
import optparse
import random
import time

# Import own general modules
import cmd_options_general
import directory

class Input_dgnm(cmd_options_general.Input,object):
    '''
    Parse arguments and store them using the optparse module.
    
     Pass sys.argv into class: Input(sys.argv).
     Options will be included in 'options' object of instance of Input_dgnm().
     Call options using:
     
     <instance>.options.<option>
     
     All other commandline arguments will be stored in <instance>.args
    '''
    
    def __init__(self,list_args, parse_parameter=True):
        # Extract scriptname to set working directory
        scriptname = list_args.pop(0)
        # If an inifile is provided, integrate values into list_args
        if ("--inifile" in list_args):
            # Insert inifile arguments before other commandline arguments.
            # This way commandline options prevail over inifile options
            print("Inifile option found in the command line.")
            list_args = self._parse_inifile(list_args) + list_args    

        usage = "usage: python %prog [options]"
        self._dbg("Argument list: %s" % list_args)
        
        # Debug print of current working directory
        self._dbg("Current workdir: %s" % os.getcwd())

        # Initialisation of the OptionParser
        parser = optparse.OptionParser(prog=scriptname,usage=usage)      
              
        # Set defaults for test mode
        parser.set_defaults(root = os.getenv("DGNM_ROOT"),
                            parameter_ini = os.path.join(os.getcwd(), "parameters.ini"),
                            species_ini = os.path.join(os.getcwd(), "species.ini"),
                            starttime = 1999,
                            endtime = 2000,
                            maskid = 99, #Rhine
                            mask_bool_operator = "EQ",
                            outputdir = os.path.join(os.getcwd(), r"../output"),
                            aggregationgrid = None,
                            water_inputdir = os.path.join(os.getcwd(), r"../hydro_netcdf"),
                            load_inputdir = os.path.join(os.getcwd(), r"../loads_netcdf"),
                            lspinup = 0,
                            allorders = 0,
                            lstrahlergrids = 0,
                            lmanning = 0,
                            lfloodplains = 0,
                            lslope = 0,
                            lbudget = 0,
                            largs = 0,
                            llocal = 1,
                            ldebug = 0
                            )
        parser.add_option("-r", "--root",
                          type = "string",
                          dest = "root",
                          action="store",
                          help="Not used anymore.")
        parser.add_option("-m", "--file_mask",
                          type = "string",
                          dest = "file_mask",
                          action="store",
                          help="Full path to ascii grid which is used as mask.")
        parser.add_option("-i", "--parameter_ini",
                          type = "string",
                          dest = "parameter_ini",
                          action="store",
                          help="Full path to file with all 'constant' parameter values like filenames etc..")                            
        parser.add_option("-s", "--species_ini",
                          type = "string",
                          dest = "species_ini",
                          action="store",
                          help="Full path to file with all definition of the species.")                            
        parser.add_option("--starttime",
                          type = "float",
                          dest = "starttime",
                          action="store",
                          help="Start year of the simulation.")
        parser.add_option("--endtime",
                          type = "float",
                          dest = "endtime",
                          action="store",
                          help="End year of the simulation.")
        parser.add_option("--maskid",
                          type = "int",
                          dest = "maskid",
                          action="store",
                          help="ID of simulated river basin used for the mask.")
        parser.add_option("--mask_bool_operator",
                          type = "string",
                          dest = "mask_bool_operator",
                          action="store",
                          help="Logical evaluation for riverid: GT, GE, LT, LE, EQ, NQ to set the mask.")
        parser.add_option("-w", "--water_inputdir",
                          type = "string",
                          dest = "water_inputdir",
                          action="store",
                          help="Input directory with the water information.")
        parser.add_option("-l", "--load_inputdir",
                          type = "string",
                          dest = "load_inputdir",
                          action="store",
                          help="Input directory with the load information for different sources.")
        parser.add_option("--fi", "--fix_input",
                          type = "string",
                          dest = "fix_input",
                          action="store",
                          help="Input directory with temporally fixed data")
        parser.add_option("-o", "--outputdir",
                          type = "string",
                          dest = "outputdir",
                          action="store",
                          help="Output directory.")
        parser.add_option("-a", "--aggregationgrid",
                          type = "string",
                          dest = "aggregationgrid",
                          action="store",
                          help="File name of grids to aggregate tabular output information on, multiple files separated by ,.")
        parser.add_option("-d","--ldebug",
                          dest = "ldebug",
                          action="store",
                          help="Print debug messages 1: yes, 0: no.")
        parser.add_option("--lspinup",
                          type = "int",
                          dest = "lspinup",
                          action="store",
                          help="Perform a spinup period at starting year. 1: yes, 0: no, -1: initialization with given species amounts.")
        parser.add_option("--lstrahlergrids",
                          type = "int",
                          dest = "lstrahlergrids",
                          action="store",
                          help="Characteristics of strahler orders are gridded. 1: yes, 0: no, use constant, e.g. Wollheim method.")
        parser.add_option("--lmanning",
                          type = "int",
                          dest = "lmanning",
                          action="store",
                          help="Calculate velocity in small streams with Manning's formula. 1: yes, 0:no-Wollheim parameters are used.")
        parser.add_option("--lfloodplains",
                          type = "int",
                          dest = "lfloodplains",
                          action="store",
                          help="Include PCR-GLOBWB floodplains as water bodies along mainstream. 1: yes, 0: no")
        parser.add_option("--lslope",
                          type = "int",
                          dest = "lslope",
                          action="store",
                          help="Include slope information in calculations. 1: yes, 0: no")
        parser.add_option("--allorders", 
                          type = "int",
                          dest = "allorders",
                          action="store",
                          help="Print outputs for lower orders. 1: yes, 0: no.")
        parser.add_option("--lbudget", 
                          type = "int",
                          dest = "lbudget",
                          action="store",
                          help="Print budget variables for all orders. 1: yes, 0: no.") 
        parser.add_option("--largs",
                          type = "int",
                          dest = "largs",
                          action="store",
                          help="Print arguments used for reactions at all timesteps. 1: yes, 0: no.") 
        parser.add_option("--llocal",
                          type = "int",
                          dest = "llocal",
                          action="store",
                          help="Perform the run on a local computer or on the cluster. 1: yes, 0: no.")                          
        parser.add_option("--inifile",
                          type = "string",
                          dest = "inifile",
                          action="store",
                          help="""Path to ini-file.
                          In inifile use commandline options, followed by '=', followed by argument, e.g. --outputdir = OECD.
                          Comments must be preceded by '#'. Everything after '#' on a line will be ignored by the parser."""
                          )

        (self.options, self.args) = parser.parse_args(args=list_args)
        # Make the same object as self.options, but with only the variable input options.
        # This is needed for writing in the log file.
        (self.options_var, self.args1) = parser.parse_args(args=list_args)


               
        # Set ldebug and lss_grw as int parameters:
        self.options.ldebug = int(self.options.ldebug)

        # Override the root parameter with the value of the environment parameter
        self.options.root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..', os.path.dirname(self.options.inifile)))
	
        # Expand all file names to an absolute reference
        self.options.species_ini = os.path.join(self.options.root, self.options.species_ini)
        self.options.parameter_ini = os.path.join(self.options.root, self.options.parameter_ini)
        self.options.outputdir = os.path.join(self.options.root, self.options.outputdir)
        
        self.options.water_inputdir = os.path.join(self.options.root, self.options.water_inputdir)
        self.options.load_inputdir = os.path.join(self.options.root, self.options.load_inputdir)

        if (self.options.aggregationgrid != None):
            self.options.aggregationgrid_number = 0
            # Make a number of attributes
            filenames = self.options.aggregationgrid.split(",")
            for item in range(len(filenames)):
                setattr(self.options,"aggregationgrid"+str(self.options.aggregationgrid_number),\
                          os.path.join(self.options.root, filenames[item]))
                self.options.aggregationgrid_number += 1                
            
        # Read all the default parameter settings out of the parameter.ini
        # Here all the other semi constant parameters for the model are set!
        self.options.parameter_ini = os.path.join(self.options.root, self.options.parameter_ini)
        self._parse_parameter_inifile(self.options,self.options.parameter_ini)
    
        # Set temp output directory for this river basin.
        # First check the temp directory
        self.options.tmp_dir = directory.ensure(os.path.join(self.options.outputdir, str(int(random.random()*1000))))     
    
        # This is to use all the directory information, which is given by the user
        if self.options.file_mask == None:
            self.options.lmask = False
        else:   
            self.options.lmask = True

        #if (self.options.lmask):
            #self.options.scenariodir = self.options.root 
            #self.options.file_mask = os.path.join(self.options.scenariodir,self.options.file_mask)

        self.options.file_mask         = os.path.join(self.options.root,self.options.file_mask)  

        # Add all the directories to the filenames
        try:
            self.options.startup_file   = os.path.join(self.options.root ,self.options.startup_file)
        except(AttributeError):
            print("No startup file is found in the parameter ini file. No start up file is used")
        self.options.ldd           = os.path.join(self.options.water_inputdir,self.options.ldd)
        self.options.basin         = os.path.join(self.options.water_inputdir,self.options.basin)
        self.options.channel_depth   = os.path.join(self.options.water_inputdir,self.options.channel_depth)
        self.options.channel_width   = os.path.join(self.options.water_inputdir,self.options.channel_width)
        self.options.slope = os.path.join(self.options.water_inputdir,self.options.slope)      
                
        # Switch debugging on/off
        self.set_debug(self.options.ldebug)
        # If no arguments are provided, print usage
        if len(list_args) == 0:
            print(parser.print_help())
            raise parser.error("No arguments provided.")
                       
        # Display all options in case of debug.
        self._dbg("Start checking the command line arguments")
        # Write all settings to screen
        self._dbg("Command line arguments found:")
        self._dbg(self.options)

        # Check given input parameters
        self.run_checks()
       

    def run_checks(self):
        '''
        Check options
        
        '''
        self.validate_directory(self.options.root, bool_write=False)
       
        # Create new output directory if necessary
        if not os.path.exists(self.options.outputdir):
            os.makedirs(self.options.outputdir)
            
        self.validate_directory(self.options.outputdir, bool_write=True)
        #self.validate_directory(self.options.scenariodir, bool_write=False)
        self.validate_directory(self.options.water_inputdir, bool_write=False)
        self.validate_directory(self.options.load_inputdir, bool_write=False)
               
        # Check the existence of all given input files                      
        if (self.options.lmask):
            self.validate_file(self.options.file_mask)
            
        self.validate_file(self.options.ldd)

        self.validate_file(self.options.basin)
        self.validate_file(self.options.channel_depth)
        self.validate_file(self.options.channel_width)
        
        # Stop checking
        return
        


