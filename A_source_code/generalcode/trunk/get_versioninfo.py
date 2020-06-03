# ******************************************************
## Revision "$LastChangedDate: 2018-07-13 09:41:57 +0200 (vr, 13 jul 2018) $"
## Date "$LastChangedRevision: 6 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/get_versioninfo.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import own general modules
from error import MyError
import get_file_info
import my_sys
import os

def print_versioninfo_datafiles(srcdir,fp=None):
    '''
    Print the information of the newest file in the srcdir.
    If pointer to file (argument fp) is equal to None then info is written to screen.
    Else it is written to the file pointer.
    '''
    file_ext = [".ASC", ".GZ", ".MAP", ".NC"]

    # Make a list of all files.
    filelist = os.listdir(srcdir)

    # List with all info
    output = []

    # Get all info from all files
    for filename in filelist:
        # Check only the files which have an extension mentioned in file_ext.
        if (os.path.splitext(filename)[1].upper() in file_ext):
            output.append(get_file_info.get_file_date_size(os.path.join(srcdir,filename)))

    # Create the info strings for this directory
    outstring = []
    if (len(output) > 0):
        outstring.append("******************************************************")

        # Sort the file information.
        output.sort(reverse=True)

        outstring.append(" Version information for directory " + srcdir)
        outstring.append(" Newest file: " + str(output[0][3]))
        outstring.append(" Last changed time: " + str(output[0][2]))
        outstring.append(" File size (in bytes): " + str(output[0][1]))

        outstring.append("******************************************************")
    else:
        outstring.append("******************************************************")
        outstring.append(" No version information found for directory " + srcdir)
        outstring.append("******************************************************")      

    # Write to screen or to file
    for line in outstring:
        if (fp == None):
            print(line)
        else:
            fp.write(line + "\n")


def get_versioninfo_files(srcdir):
    '''
    Read the header of all python scripts and returns the highest number of any file as results.
    The highest number is determined by LastChangedRevision.
    Version number must be in the first 100 lines.
    '''

    versionnumber = 0
    filename_versionnumber = None

    # Make a list of all files.
    filelist = os.listdir(srcdir)

    for filename in filelist:
        # Check only the python files.
        if (os.path.splitext(filename)[1].upper() == ".PY"):
            lines = my_sys.my_readfile(os.path.join(srcdir,filename))
            for item in range(len(lines)):
                if (item > 100):
                    break
                if ("$LastChangedRevision:" in lines[item]):
                    # Get version number. Information is given between the $ signs.
                    ind_end = lines[item].rfind("$")
                    ind_start = lines[item].find("$")
                    version = int(lines[item][ind_start+1:ind_end].split(":")[1])
                    if (version > versionnumber):
                        versionnumber = version
                        filename_versionnumber = filename
                    break

    return versionnumber, filename_versionnumber

def print_versioninfo(filename,fp=None):
    '''
    if pointer to file (argument fp) is equal to None then info is written to screen.
    Else it is written to the file pointer.
    Version info is fetched from filename.
    '''
    labels = ["$LastChangedDate:","$LastChangedRevision:","$LastChangedBy:","$HeadURL:"]
    lfound = len(labels) * [False]

    # Read filename
    lines = my_sys.my_readfile(filename)

    # Open version information
    if (fp == None):
        print("******************************************************")
        print(("## Version information from file: " + filename))
    else:
        fp.write("******************************************************")
        fp.write("## Version information from file: " + filename)

    # Get version information from file
    for line in lines:
        for item in range(len(labels)):
            if (not lfound[item]):
                if (labels[item] in line):
                   if (fp == None):
                       print(line)
                   else:
                       fp.write(line)
                   lfound[item] = True
                   break

    # Close version information
    if (fp == None):
        print("******************************************************")
    else:
        fp.write("******************************************************")

def get_versioninfo(obj,outputdir):
    '''
    Get version (svn) information of directories and return a list of strings.
    '''
    try:
        my_sys.my_system(os.path.join("../../tools/sliksvn/bin","svn") + " info " + str(obj) + " 1> " +\
                         os.path.join(outputdir,"tempuit1.txt") + " 2> " +\
                         os.path.join(outputdir,"tempuit.err"),fatal=1)
        my_sys.my_system(os.path.join("../../tools","grep") + " Rev " +\
                         os.path.join(outputdir,"tempuit1.txt") + " 1> " +\
                         os.path.join(outputdir,"tempuit.txt") + " 2> " +\
                         os.path.join(outputdir,"tempuit.err"),fatal=1)
        list = my_sys.my_readfile(os.path.join(outputdir,"tempuit.txt"))
        my_sys.my_removefile(os.path.join(outputdir,"tempuit.txt"))
        my_sys.my_removefile(os.path.join(outputdir,"tempuit1.txt"))
        my_sys.my_removefile(os.path.join(outputdir,"tempuit.err"))
        list.insert(0,"Directory " +str(obj))
    except (MyError,OSError):
        try:
            #Try another time but with the assumption that the system has svn installed.
            my_sys.my_system("svn info " + str(obj) + " 1> " +\
                             os.path.join(outputdir,"tempuit1.txt") + " 2> " +\
                             os.path.join(outputdir,"tempuit.err"),fatal=1)
            my_sys.my_system("grep Rev " +\
                             os.path.join(outputdir,"tempuit1.txt") + " 1> " +\
                             os.path.join(outputdir,"tempuit.txt") + " 2> " +\
                             os.path.join(outputdir,"tempuit.err"),fatal=1)
            list = my_sys.my_readfile(os.path.join(outputdir,"tempuit.txt"))
            my_sys.my_removefile(os.path.join(outputdir,"tempuit.txt"))
            my_sys.my_removefile(os.path.join(outputdir,"tempuit1.txt"))
            my_sys.my_removefile(os.path.join(outputdir,"tempuit.err"))
            list.insert(0,"Directory " +str(obj))
        except (MyError,OSError): 
            list = []
            list.append("Directory "+str(obj) + " has no version info.")
            # Remove temporary files which are created here 
            if (os.path.isfile(os.path.join(outputdir,"tempuit.txt"))):
                my_sys.my_removefile(os.path.join(outputdir,"tempuit.txt"))
            if (os.path.isfile(os.path.join(outputdir,"tempuit1.txt"))):
                my_sys.my_removefile(os.path.join(outputdir,"tempuit1.txt"))
            if (os.path.isfile(os.path.join(outputdir,"tempuit.err"))):
                my_sys.my_removefile(os.path.join(outputdir,"tempuit.err"))

    return list
