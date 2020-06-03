# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/test1.py $"
# ******************************************************
import os
import sys
import random

def plus(array,item):
    tel = random.randint(0,10000)
    dum = tel
    for j in range(tel):
        dum += 1
    return array[item]

def plusarray(array,ib,ie):
    uit = []
    for j in range(ib,ie):
        uit.append(plus(array,j))
    return uit

if __name__ == "__main__":

    import pp
    job_server = pp.Server()
    start = 0
    endprog = 500
    parts = 50
    step = (endprog - start) / parts + 1

    job_server.set_ncpus()
    print(job_server.get_ncpus())
    array=[]
    for item in range(start,endprog):
        array.append(item)
    jobs=[]
    out=[]
    for index in range(parts):
        starti = start+index*step
        endi = min(start+(index+1)*step, endprog)
        print(starti,endi)
        jobs.append(job_server.submit(plusarray,(array,starti,endi),(plus,),("random",)))
        #for job in jobs:
            #print job()
            #out.extend(job())
        #out = [job() for job in jobs]
        #print job_server.print_stats()
        #plusarray(array,starti,endi)
    #print out
    for job in jobs:
        out.extend(job())
    print(out)
