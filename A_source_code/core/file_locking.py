# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/file_locking.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os
import time
import random
import sys
import traceback

def aquire_lock(params,lock,procid,locktext=None):
    '''
    Aquire the lock, or when the lock is None, then 
    use own locking process to aquire the lock (own made).
    With the locktext you can  have more than one lock. 
    '''
    if (lock != None):
        # Here the locking of python is used.
        lock.acquire()
    else:
        # The locking file is stored on the output directory.
        if (locktext == None):
            filename = os.path.join(params.outputdir,"lock_general.lock")
            filename1 = os.path.join(params.outputdir,"lock_general1.lock")
        else:
            filename = os.path.join(params.outputdir,"lock_"+str(locktext)+".lock")
            filename1 = os.path.join(params.outputdir,"lock_"+str(locktext)+"1.lock")

        laquire = False
        counter = 0
        lfirst = True
        while not (laquire):
            lowner = False
            try:
                # Check whether locking file already exists, which means that there is another process already active.
                if (not os.path.exists(filename1)):
                    try:
                        # First create this file.
                        with open(filename1, "w") as fp:
                            fp.write(str(0))
                        fp.close()
                        lowner = True
                    except:
                        print("*2a*****************************************************************")
                        print(sys.exc_info()[0])
                        print(traceback.print_exc())
                        time.sleep(0.3)
                        #pass                        

                    try:
                        if (lowner):
                            # Wait small amount of time to make sure that there are not other processes that also have been writen to this lock file.
                            time.sleep(random.uniform(0.3,0.6))
                            print("Sleep max of 0.6 sec") 

                            # First write our proces-id into this file.
                            if (not os.path.exists(filename)):
                                with open(filename, "w") as fp:
                                    fp.write(str(procid))
                                fp.close()

                                # Wait small amount of time to make sure that there are not other processes that also have been writen to this lock file.
                                time.sleep(0.4)
                                print("Sleep of 0.4 sec") 

                                # Check whether this file is the file of this process.
                                try:
                                    with open(filename, "r") as fp:
                                        if (str(procid) == fp.readline()):
                                            laquire = True
                                    fp.close()
                                except:
                                    # Other process could have thrown away the file
                                    print("*1*****************************************************************")
                                    print(sys.exc_info()[0])
                                    print(traceback.print_exc())
                                    #pass
                    except:
                        print("*2*****************************************************************")
                        print(sys.exc_info()[0])
                        print(traceback.print_exc())
                        time.sleep(0.2)
                        #pass

                elif (counter == 1):
                    if (os.path.exists(filename)):
                        # Check whether this file is the file of this process.
                        try:
                            with open(filename, "r") as fp:
                                if (str(procid) == fp.readline()):
                                    laquire = True
                            fp.close()
                        except:
                            # Other process could have thrown away the file
                            print("*3*****************************************************************")
                            print(sys.exc_info()[0])
                            print(traceback.print_exc())                            
                            #pass
                        if (not laquire):
                            time.sleep(0.37)
                            print("Sleep of 0.37")
                    else:
                        counter = 0
                else:
                    # File already exists, wait 0.4 seconds.
                    time.sleep(0.4)
                    if (lfirst):
                        print("Sleep of 0.4")
                        lfirst = False
            except:
                # Now the job knows that there is something wrong.
                print("*********** Something went wrong with file_locking. *****************************************")
                print("*4*****************************************************************")
                print(sys.exc_info()[0])
                print(traceback.print_exc())
                time.sleep(0.2)
                counter = 1


def release_lock(params,lock,locktext=None):
    '''
    Release the lock, or when the lock is None, then 
    use own locking process to release the lock (own made).
    With the locktext you can  have more than one lock. 
    '''
    if (lock != None):
        # Here the locking of python is used.
        lock.release()
    else:
        # The locking file is stored on the output directory.
        if (locktext == None):
            filename = os.path.join(params.outputdir,"lock_general.lock")
            filename1 = os.path.join(params.outputdir,"lock_general1.lock")
        else:
            filename = os.path.join(params.outputdir,"lock_"+str(locktext)+".lock")
            filename1 = os.path.join(params.outputdir,"lock_"+str(locktext)+"1.lock")

        # Throw away the locking file.
        if (os.path.exists(filename)):
            try:
                os.remove(filename)
                os.remove(filename1)
            except:
                print("*7*****************************************************************")
                print(sys.exc_info()[0])
                print(traceback.print_exc())             
                while (os.path.exists(filename)):
                    try:
                        os.remove(filename)
                    except:
                        print("Try again to remove file " + filename + " with small sleep 0.05")
                        print("*5*****************************************************************")
                        print(sys.exc_info()[0])
                        print(traceback.print_exc())                            
                        time.sleep(0.05)
                while (os.path.exists(filename1)):
                    try:
                        os.remove(filename1)
                    except:
                        print("Try again to remove file " + filename1 + " with small sleep 0.05")
                        print("*********** Something went wrong with file_locking. *****************************************")
                        print("*6*****************************************************************")
                        print(sys.exc_info()[0])
                        print(traceback.print_exc())     
                        time.sleep(0.05)
   


def clean_locks(params):
    '''
    Clean the directory of all locks that are placed on the output directory.
    This is a routine needed before the start of the processes.
    '''
    list = os.listdir(params.outputdir)
    for ifile in list:
        if (os.path.isfile(os.path.join(params.outputdir,ifile))):
            if ifile.endswith(".lock"):
               os.remove(os.path.join(params.outputdir,ifile))

