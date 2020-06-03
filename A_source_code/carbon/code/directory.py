# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import glob
import os

def get_files(searchdir):
    '''
    returns all files in directory "searchdir" as a list
    '''
    return [f for f in os.listdir(searchdir) if os.path.isfile(os.path.join(searchdir, f))]


def get_files_of_type(searchdir, ext):
    '''
    returns all files of type *.ext in a directory "searchdir" as a list
    '''
    ext = "."+ext.replace('.', '')
    file_list = glob.glob(os.path.join(searchdir, '*'+ ext))
    return file_list

def get_files_with_str(searchdir, searchstr, exclude=[]):
    '''
    returns all files containing string "searchstr" in directory "searchdir"
    and excludes files containing strings in list exclude

    wildcards are advised
    '''
    file_list = list()
    searchdir=searchdir+'/'
    for filename in glob.iglob(os.path.join(searchdir,searchstr), recursive=True):
      file_list.append(filename)
    exclude_list = []
    for str in exclude:
      for file in file_list:
        if str in file:
          exclude_list.append(file)
    for item in exclude_list:
      try:
        file_list.remove(item)
      except(ValueError):
        pass
      
    return file_list
  
def get_dirs(searchdir, full=False):
    '''
    returns all directory names in directory "searchdir"

    full = False returns only the name of the directory
    full = True returns the entire path of directory
    '''  
    if full==True:
      dirs = list()  
      dirs_dum = [x[1] for x in os.walk(searchdir)][0]
      for d in dirs_dum:
          dirs.append(os.path.join(searchdir, d))
    else:
      dirs = [x[1] for x in os.walk(searchdir)][0]  
    return dirs
  
def ensure(d, returnB=True):
    '''
    function to make sure directory "d" exists

    returnB = False makes the function return nothing
    returnB = True makes the function return the ensured directory "d"
    '''
    if not os.path.exists(d):
        os.makedirs(d)
    if returnB:
        return d   

def up(d):
    '''
    returns one directory above directory "d" as entire path
    '''
    return os.path.abspath(os.path.join(d, os.pardir))

def get_parent(d, full=False):
    '''
    returns one directory above directory "d"

    full = False returns only the name of the returned directory
    full = True returns the entire path of returned directory
    '''
    if full:
        return up(d)
    else:
        return os.path.basename(os.path.normpath(d))