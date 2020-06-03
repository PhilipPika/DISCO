# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/make_cluster_cmd.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import general python modules
import os
import pickle

# Import own general modules
from error import *
import get_versioninfo
import make_slurm_script
import my_logging
import my_sys

# Import local modules.
import directory
import general_func

def make_cluster_cmd(params,jobnum,numbers,total_cells):
    '''
    Create a string which can be used to start the job on a cluster.
    '''        
    # Dump river numbers to file.
    fp_out = open(os.path.join(params.outputdir,"rivers_"+str(jobnum)+".pkl"),"wb")
    pickle.dump(numbers,fp_out, -1)
    fp_out.close()
            
    # Check whether there is a need for conversion of the paths, because this is used in a dos prompt or linux terminal
    outputdir = general_func.path_conversion(params.outputdir)

    # Make extra text for the submission of the jobs
    # Assign more cpu's to the job
    ncores = general_func.get_number_of_cores(params,total_cells)

    # Nodename
    try:
        nodestring = params.nodes
    except AttributeError:
        nodestring = "1"

    # Set core and node text.
    if (params.submission_text.startswith("bsub")):
        coretxt = " -n" + str(ncores) + " "
        if (nodestring != "1"):
            coretxt += " -m " + nodestring + " "
    elif (params.submission_text.startswith("qsub")):
        coretxt = " -l nodes=" + nodestring + ":ppn=" + str(ncores) + " "
    elif (params.submission_text.startswith("sbatch")):
        nodes=nodestring
        tasks=str(ncores*2)
        tasktxt = " --ntasks="+tasks+ " "
        nodetxt = " --nodes="+nodes+ " "
            
    # ERROR file
    errtxt = " -e " + os.path.join(outputdir,"error_"+str(jobnum)+".txt") + " "
    # STDOUT file
    stdtxt = " -o " + os.path.join(outputdir,"stdout_"+str(jobnum)+".txt") + " "

    # JOBname

    if (params.submission_text.startswith("bsub")):
        jobtxt = " -J "
    elif (params.submission_text.startswith("qsub")):
        jobtxt = " -N "
    elif (params.submission_text.startswith("sbatch")):
        jobtxt = " --job-name "
    jobtxt += os.path.splitext(params.inifile)[0] +"_" +str(jobnum) + "_"+"dgnm " 
            
    # Put extra text after bsub command.
    fields = params.submission_text.split()
    if (params.submission_text.startswith("bsub")):
        submission_text = fields[0] + errtxt + jobtxt + coretxt + stdtxt + " ".join(fields[1:])
    elif (params.submission_text.startswith("qsub")):
        submission_text = params.submission_text + errtxt + jobtxt + coretxt + stdtxt + " -d " + os.curdir

    # Make path to main script
    if "dgnm_river.py" in directory.get_files(os.path.join(os.getenv("DGNM_ROOT"),"users", os.getenv("DGNM_USER"), "code")):
      mainscript = os.path.join(os.getenv("DGNM_ROOT"),"users", os.getenv("DGNM_USER"), "code", "dgnm_river.py")
    else:
      mainscript = os.path.join(os.getenv("DGNM_ROOT"),"core","dgnm_river.py")

    # Make the basic command.
    if (params.submission_text.startswith("bsub")):
        if (not params.submission_text.endswith("python")):
            submission_text += " python "
        cmd = submission_text + " " + mainscript + " " +\
              os.path.join(outputdir,"general_objects.pkl") + " " +\
              os.path.join(outputdir,"rivers_"+str(jobnum)+".pkl")
    elif (params.submission_text.startswith("qsub")):
        cmd = "echo python " + mainscript + " " +\
               os.path.join(outputdir,"general_objects.pkl") + " " +\
               os.path.join(outputdir,"rivers_"+str(jobnum)+".pkl") + " " +\
               submission_text
    elif (params.submission_text.startswith("sbatch")):
        cmd_str = "python " + mainscript + " " + os.path.join(outputdir,"general_objects.pkl") + " " + os.path.join(outputdir,"rivers_"+str(jobnum)+".pkl")
        # make sbatch.sh for the model run
        make_slurm_script.do(cmd_str, tasktxt, nodetxt, errtxt, stdtxt, jobtxt, outputdir)
        cmd = params.submission_text+' '+os.path.join(outputdir, 'sbatch.sh')
        #print(cmd)
    return cmd
