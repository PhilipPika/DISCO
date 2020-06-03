# ******************************************************
## Revision "$LastChangedDate: 2018-07-11 16:25:27 +0200 (wo, 11 jul 2018) $"
## Date "$LastChangedRevision: 4 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/get_file_info.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os, time
from stat import * # ST_SIZE etc

def get_file_info(filename):
    '''
    Gets information of a file and puts it in a list of strings.
    '''
    out = []
    try:
       st = os.stat(filename)
    except IOError:
        print("** WARNING ** Failed to get information about " + filename)
    else:
        out.append("********************************************")
        out.append("File name: " + os.path.abspath(filename))
        out.append("File size: " + str(st[ST_SIZE]))
        out.append("File modified: " + str(time.asctime(time.localtime(st[ST_MTIME]))))
        out.append(" ")
    return out

def get_file_date_size(filename):
    '''
    Returns date and size information of a file.
    '''
    try:
       st = os.stat(filename)
    except IOError:
        print("** WARNING ** Failed to get information about " + filename)
        file_time = 0
        file_size = None
        file_date = None
    else:
        file_time = st[ST_MTIME]
        file_size = st[ST_SIZE]
        file_date = time.asctime(time.localtime(st[ST_MTIME]))
    
    return file_time, file_size, file_date, filename
