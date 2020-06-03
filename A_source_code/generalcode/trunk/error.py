# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/error.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
## ---------------------------------------------------------------------------
## -
## -    Name:    error.py
## -    Author:  Arthur Beusen en Martine de Vos PBL/IMP
## -    Adjusted:Allard Warrink PBL/IMP
## -    Date:    June 18 2008
## -    Purpose: Error handling of exceptions
## -
## ---------------------------------------------------------------------------

class MyError(Exception):
    """
    Stores an exception on type or format of arguments from the commandline. 
    This class is used to generate a message to standard output.
    This message can contain multiple lines.
    Inherits properties from build in class Exception.
    """
    def __init__(self, *args):
        """ Initialisation of class MyError."""
        self.args = args
    
    def __str__(self):
        """ Returns a formatted version of the objects content. """
        return repr(''.join(map(str, self.args)))
    
    def write(self):
        """ Writes an error message to standard output. """
        print("***** ERROR ******")
        for arg in self.args:
            print(arg)

class MyError0(MyError):
    '''
    For backward compatibility
    '''
    def write(self):
        """ Writes an error message to standard output. """
        for arg in self.args:
            print(arg)


class MyError1(MyError):
    '''
    For backward compatibility
    '''
    pass

class MyError2(MyError):
    '''
    For backward compatibility
    '''
    pass
