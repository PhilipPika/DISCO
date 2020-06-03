# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/maintest.py $"
# ******************************************************
# This code is taken from: http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method#Python  on 15 november 2011

import os
import sys


if __name__ == "__main__":

    # Import Global ("generalcode") modules
    # Generic code files ideally should be placed in the directory "generalcode"
    # subfolder, but they may also be placed at the same base folder
    __general = os.path.join(r'..')
    if os.path.exists(__general): 
        sys.path.insert(0, __general)
        print(__general + " is added to the python search path for modules.") 
    import arthur
    import ppsimpson
    #from math import sin,cos,log,exp
    from math import *
    import time
    import pp
    job_server = pp.Server()
    start = 1
    endprog = 50
    parts = job_server.get_ncpus()
    step = (endprog - start) / parts + 1

    job_server.set_ncpus()
    print(job_server.get_ncpus())

    print("Estim solution: ",ppsimpson.adaptive_simpsons_rule(sin,0,1,.000000001))
    print("Exact solution: ",-cos(1) + cos(0)) 
    print("Estim solution: ",ppsimpson.adaptive_simpsons_rule(sin,0,1,.001)) 
    print("Exact solution: ",-cos(1) + cos(0)) 
    for ncpu in (2,16):
        start_time=time.time()
        job_server.set_ncpus(ncpu)
        jobs=[]
        for index in range(parts):
            starti = start+index*step
            endi = min(start+(index+1)*step, endprog)
            print(starti,endi)
            #jobs.append(job_server.submit(adaptive_simpsons_rule_array,(arthur.arthur,0,1,0.0000000000001,starti,endi),(recursive_asr,simpsons_rule))))
            jobs.append(job_server.submit(ppsimpson.adaptive_simpsons_rule_array,(arthur.arthur,0,1,0.0000000000001,starti,endi),(ppsimpson.adaptive_simpsons_rule,ppsimpson.recursive_asr,ppsimpson.simpsons_rule),("arthur","ppsimpson")))
            #jobs.append(job_server.submit(sin,(.1)))
        out = [job() for job in jobs]
        #print out
        end_time = time.time()
        print("Time elapsed for " +str(ncpu) + " workers in  " +str(end_time-start_time) + "sec")
        print(job_server.print_stats())
