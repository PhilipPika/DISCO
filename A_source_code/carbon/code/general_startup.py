# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import general python modules
import os

# Import own general modules
from error import *
import get_versioninfo
import my_logging
import my_sys

# Import local modules.
import cmd_options_dgnm

def general_startup(args):

    # Parse command-line arguments and set parameters for script
    try:
        param = cmd_options_dgnm.Input_dgnm(args)
        params = param.options
        params_var = param.options_var
    except SystemExit:
        raise MyError("Error has occured in the reading of the commandline options.")  
    
    # Start timer and logging
    s = my_sys.SimpleTimer()
    log = my_logging.Log(params.outputdir,"%i_%i.log" % (params.starttime,params.endtime))
    print("Log will be written to %s" % log.logFile)
    
    # If no arguments are provided do a run with defaults in params
    if len(args) == 0:
        log.write_and_print("No arguments provided: starting default run...")
        
    # Write parameter settings and version information to log file
    if params.llocal==1:
      log.write_and_print("Starting run....")
    else:
      log.write_and_print("Starting submission....")
    log.write("#############################################",print_time=False,lcomment=False)
    log.write("# Ini input options used:",print_time=False,lcomment=False)
    for option in str(params_var).split(", "):
        log.write(" --%s = %s" % (option.split(": ")[0].strip("{").strip("'"),
                               option.split(": ")[1].strip("}").strip("'")),
                print_time=False,lcomment=False)
    log.write("#############################################",print_time=False,lcomment=False)               
    log.write("# All input options used:",print_time=False,lcomment=False) 
    for option in str(params).split(", "):
        log.write("# %s = %s" % (str(option).split(": ")[0].strip("{").strip("'"),
                               str(option).split(": ")[1].strip("}").strip("'")),
                print_time=False,lcomment=False)
    log.write("#############################################",print_time=False,lcomment=False)                
      
    # Check whether there are command-line arguments which are not used
    if (len(param.args) > 0):
        txt = "The following command line arguments will not be used:"
        log.write_and_print(txt + str(param.args))

    # Write svn information of input and scripts to log file.

    log.write("******************************************************",print_time=True,lcomment=True)
    dirname = os.path.join(os.getenv("DGNM_ROOT"),"core")    
    log.write("# Version information of main script.",print_time=True,lcomment=True)
    get_versioninfo.print_versioninfo(os.path.join(dirname,"dgnm_main.py"),log)

    # Write version information of generalcode directory
    dirname = os.getenv("DGNM_GENERALCODE")
    versionnumber, filename_versionnumber = get_versioninfo.get_versioninfo_files(dirname)
    log.write("# Version information of generalcode directory.",print_time=True,lcomment=True)
    get_versioninfo.print_versioninfo(os.path.join(dirname,filename_versionnumber),log)

    # Write version information of corecode directory
    dirname = os.path.join(os.getenv("DGNM_ROOT"),"core")
    versionnumber, filename_versionnumber = get_versioninfo.get_versioninfo_files(dirname)
    log.write("# Version information of core code directory.",print_time=True,lcomment=True)
    get_versioninfo.print_versioninfo(os.path.join(dirname,filename_versionnumber),log)

    # Write version information of corecode directory
    dirname = os.path.join(os.getenv("DGNM_ROOT"),"users",os.getenv("DGNM_USER"),"code")
    if (os.path.isdir(dirname)):
        versionnumber, filename_versionnumber = get_versioninfo.get_versioninfo_files(dirname)
        log.write("# Version information of personal code directory.",print_time=True,lcomment=True)
        get_versioninfo.print_versioninfo(os.path.join(dirname,filename_versionnumber),log)
    else:
        log.write("# No version information of personal code directory.",print_time=True,lcomment=True)


    return params,log,s
