# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/my_logging.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os
import time
import getpass

class Log(object):
    """Logfile class"""
    def __init__(self,logdir,filename):
        """Create and open logfile"""
        #fileName = "Log_%s_" % scenario + time.strftime("%Y%m%dT%H.%M.%S") + ".txt"
        self.logFile = os.path.join(logdir,filename)
        self.log = open(self.logFile, 'w')
        self.write("Logfile started by user: %s" % getpass.getuser())

    def write(self,message,print_time=True,lcomment=True):
        """Write message to logfile and start new line."""
        if print_time and lcomment:
            self.log.write("# %-26s: %s\n" % (time.asctime(),message))
        elif print_time and not lcomment:
            self.log.write("%-26s: %s\n" % (time.asctime(),message))
        elif not print_time and lcomment:
            self.log.write("# %-26s: %s\n" % ("",message))
        else:
            self.log.write("%s\n" % (message))
    
    def write_and_print(self,message,print_time=True):
        if print_time:
            print("# %-26s: %s" % (time.asctime(),message))
            self.write(message,print_time)
        else:
            print("# %-26s: %s" % ("",message))
            self.write(message,print_time)
    
    def __del__(self):
        self.log.close()
