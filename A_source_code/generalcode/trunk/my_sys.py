# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 16:48:30 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 2 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/my_sys.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

##Bestandsnaam:     my_sys.py
##Gemaakt:          08-11-2007
##Auteurs:          Arthur Beusen en Martine de Vos, Arjan van der Put
##Beschrijving:     Generieke module met functies voor standaard operaties op files en directories
##Functies:         my_rmdir,my_copyfile,my_copytree,my_copydir,my_removefile,
##                  my_system,my_readfile,my_getcolumn,my_getcolumnnumbers
##Classes:          SimpleTimer


##NB1: Source code om de verschillende functies te testen staat onderaan het script
##NB2: Debug informatie wordt geprint als my_sys_ldebug = 1, dit kan als argument 
##     met de functies worden meegegeven. Default waarde voor my_sys_ldebug is 0
##NB3: De fatal-optie geeft de gebruiker de mogelijkheid om in het hoofdprogramma een foutmelding te
##     genereren of zelfs het hoofdprogramma te beeindigen)mocht de geimporteerde functie fout lopen.

# Dit is een test 

import os
import time
import sys
import shutil
from error import *
import subprocess

class SimpleTimer(object):
    """Record time gone by to get an indication of performance."""
    def __init__(self):
        # start timing
        self.starttime = time.clock()
        self.lasttime = self.starttime
        
    def print_interval(self, desc):
        """print time gone by since last call to function
        
        desc        string describing the action performed between calls
        
        """
        print(self.interval(desc))
        
    def interval(self, desc):
        """return string containing time gone by since last call to function
        
        desc        string describing the action performed between calls
        
        """
        interval = time.clock() - self.lasttime
        self.lasttime = time.clock()
        return "%s took %f seconds to complete." % (desc,interval)
    
    def print_total(self, desc):
        """print time gone by since initiation
        
        desc        string describing action since initiation
        
        """
        print(self.total(desc))
        
    def total(self, desc):
        """return string containing time gone by since initiation
        
        desc        string describing action since initiation
        
        """
        total_time = time.clock() - self.starttime
        return "%s took %f seconds to complete." % (desc,total_time)
        
    def reset(self):
        self.__init__()


def my_print(text,fatal=0):
    """Prints a text to screen or raise and MyError with text as message"""
    if (fatal == 1):
        raise MyError(text)
    else:
        print(text)
 
def my_rmdir(dirname,my_sys_ldebug=0,fatal=0):
    """
    Verwijdert een directory, inclusief alle bijbehorende files en subdirectories.
    Argumenten: pad naar de te verwijderen directory en evt een debug-optie en een fatal-optie.
    NB: 
    """
    try:    
        # Controleer of opgegeven pad een bestaande directory is, zo niet: stap uit functie 
        if not os.path.isdir(dirname):
            my_print("MY_RMDIR: " + dirname + " is niet een bestaande directory.",1)

        curdir = dirname
        if (my_sys_ldebug): 
            print("MY_RMDIR: CURDIR: ",curdir)
        list = os.listdir(dirname)
        for item in range(0,len(list)):    
            if os.path.isdir(os.path.join(curdir,list[item])):
                if (my_sys_ldebug):
                    print("MY_RMDIR: DIR ",os.path.join(curdir,list[item])," wordt leeg gemaakt.")
                # Als functie op een subdirectory stuit,
                # wordt functie recursief aangeroepen (met nieuwe dirname). 
                my_rmdir(os.path.join(curdir,list[item]),my_sys_ldebug)
            else:
                if (my_sys_ldebug): print("MY_RMDIR: FILE ",os.path.join(curdir,list[item])," wordt verwijderd.")
                # Files binnen (sub)directory worden verwijderd
                try:
                    os.remove(os.path.join(curdir,list[item]))
                except:
                    my_print("MY_RMDIR: FILE " + str(os.path.join(curdir,list[item])) + " kan NIET worden verwijderd.",0) 
        if os.path.isdir(dirname):
            # Lege(sub)directory wordt verwijderd
            try:
                os.rmdir(dirname)
            except:
                my_print("MY_RMDIR: DIR" + str(dirname) + " kan NIET worden verwijderd.",0)
                
        if (my_sys_ldebug): print("MY_RMDIR: DIR ",dirname," wordt verwijderd.")
    except:
        my_print("MY_RMDIR: Verwijderen van DIR "+dirname+" is mislukt.",fatal)


    
def my_copyfile(source,destination,my_sys_ldebug=0,fatal=1):
    """
    Kopieert een bestand van source naar destination met behoud van file attributes.
    Argumenten: paden van source en destination evt een debug-optie en een fatal-optie.
    NB: Verwijdert altijd het oude bestand voordat de kopie wordt neergezet.
    """
    # Controleer of opgegeven source een bestaand bestand is, zo niet: stap uit functie
    try:
        if not os.path.isfile(source):
            my_print("MY_COPYFILE: " + str(source) + " is niet een bestaand bestand.",1)

        if not os.path.isdir(destination):
            if os.path.isfile(destination):
                if (my_sys_ldebug): print("MY_COPYFILE: Oude FILE ",destination," wordt verwijderd.")
                 # Als er nog een oud bestand staat, wordt dit verwijderd.
                 # Als dit mislukt: stap uit functie
                try:
                    os.remove(destination)
                except:
                    my_print("MY_COPYFILE: Oude FILE " + str(destination) + " kan NIET worden verwijderd.",1)

                # Nieuwe bestand wordt gekopieerd. Als dit mislukt: stap uit functie
                # Preserve file attributes
            try:
                shutil.copy2(source,destination)
            except:
                my_print("MY_COPYFILE: FILE "+ str(source) + " kan NIET worden gekopieerd naar " + str(destination),1)
       
            if (my_sys_ldebug):
                print("MY_COPYFILE: FILE ",destination," is gekopieerd.")
        else:
            if (my_sys_ldebug):
                print("MY_COPYFILE: FILE "+source+" wordt naar "\
                      + destination+"\\"+os.path.split(source)[1] + " gekopieerd.")
            # Als destination een directory is,
            # wordt functie recursief aangeroepen (met nieuwe destination).
            my_copyfile(source,os.path.join(destination,os.path.split(source)[1]),my_sys_ldebug)
    except:
        my_print("MY_COPYFILE: Kopieren van FILE " + str(source) + " naar " + str(destination) + " is mislukt.",fatal)


def my_copytree(source,destination,my_sys_ldebug=0,fatal=1):
    """Kopieert alleen de directory structuur van source naar destination
       Argumenten: paden van source en destination evt een debug-optie en een fataloptie.
    """
    try:
        # Controleer of opgegeven source een bestaande directory is, zo niet: stap uit functie.
        if not os.path.isdir(source):
            my_print("MY_COPYTREE: " + str(source) + " is niet een bestaande directory.",1)

        if not os.path.isdir(destination):
            if (my_sys_ldebug):
                print("MY_COPYTREE: DIR " + destination + " wordt gemaakt.")

            try:
                os.makedirs(destination)
            except:
                my_print("MY_COPYTREE: DIR " + str(destination) \
                      + " kan niet worden aangemaakt.",1)

        list = os.listdir(source)
        for item in range(0,len(list)):
            # Controleer of source directory nog een subdirectory heeft,
            # zo ja: spring naar subdirectory.  
            if (not os.path.isfile(os.path.join(source,list[item]))):
                # Controleer of subdirectory al bestaat in destination, zo niet maak subdirectory.     
                if not os.path.isdir(os.path.join(destination,list[item])):
                    if (my_sys_ldebug):
                        print("MY_COPYTREE: DIR ",os.path.join(destination,list[item])+" wordt gemaakt.")
                    # Nieuwe subdirectory wordt aangemaakt.
                    try:
                        os.makedirs(os.path.join(destination,list[item]))
                    except:
                        my_print("MY_COPYTREE: DIR" + str(os.path.join(destination,list[item])) \
                              +" kan niet worden aangemaakt.",fatal)
                
                # Vanuit subdirectory in source wordt functie recursief aangeroepen
                # (met nieuwe source en destination).     
                my_copytree(os.path.join(source,list[item]),os.path.join(destination,list[item]),my_sys_ldebug)
    except:
        my_print("MY_COPYTREE: Kopieren van tree " + str(source) + " naar " + str(destination) + " is mislukt.",fatal)

    
def my_copydir(source,destination,my_sys_ldebug=0,fatal=1):
    """
    Kopieert een directory van source naar destination, inclusief alle files en subdirectories
    Argumenten: paden van source en destination en evt een debug-optie en een fatal-optie.
    NB: Maakt de nieuwe directory structuur mbv my_copytree, kopieert files mbv my_copyfile
    """
    try:
        # Directory structuur wordt gekopieerd mbv my_copytree.
        my_copytree(source,destination,my_sys_ldebug,fatal)

        list = os.listdir(source)
        for item in range(0,len(list)):
            # Controleer of source files heeft in (sub)directories, zo ja: kopieer deze .
            if os.path.isfile(os.path.join(source,list[item])):
                if (my_sys_ldebug):
                    print("MY_COPYDIR: FILE: " + os.path.join(source,list[item])\
                          + " wordt naar my_copyfile gestuurd.")
                my_copyfile(os.path.join(source,list[item]),destination,my_sys_ldebug,fatal)
            else:
                if (my_sys_ldebug):
                    print("MY_COPYDIR: DIR: " + os.path.join(source,list[item])\
                          + " wordt naar my_copydir gestuurd.")
                # Vanuit subdirectory in source wordt functie recursief aangeroepen
                #(met nieuwe source en destination). 
                my_copydir(os.path.join(source,list[item]),os.path.join(destination,list[item]),my_sys_ldebug)
    except:
        my_print("MY_COPYDIR: Kopieren van DIR " + str(source) + " naar " + str(destination) + " is mislukt.",fatal)

def my_copydir1(source,destination,my_sys_ldebug=0,fatal=1):
    """
    Kopieert alle bestanden van een directory van source naar destination, exclusief de subdirectories
    Argumenten: paden van source en destination en evt een debug-optie en een fatal-optie.
    NB: Maakt de nieuwe directory structuur mbv my_copytree, kopieert files mbv my_copyfile
    """
    try:
        # Directory structuur wordt gekopieerd mbv my_copytree.
        list = os.listdir(source)
        for item in range(0,len(list)):
            # Controleer of source files heeft in (sub)directories, zo ja: kopieer deze .
            if os.path.isfile(os.path.join(source,list[item])):
                if (my_sys_ldebug):
                    print("MY_COPYDIR: FILE: " + os.path.join(source,list[item])\
                          + " wordt naar my_copyfile gestuurd.")
                my_copyfile(os.path.join(source,list[item]),destination,my_sys_ldebug,fatal)
    except:
        my_print("MY_COPYDIR: Kopieren van DIR " + str(source) + " naar " + str(destination) + " is mislukt.",fatal)



def my_system(cmd,my_sys_ldebug=0,fatal=1):
    """
    Doet een systemcall
    Argumenten: dos commando en evt een debug-optie en een fatal-optie
    """
    error=0
    try:    
        error=subprocess.call(cmd,shell=True)
        if (my_sys_ldebug):
            print("MY_SYSTEM: System call: "+ cmd + " has finished.")
        if (error < 0):
            my_print("MY_SYSTEM: System call is terminated met signal: " + str(error),0)
            my_print("MY_SYSTEM: System call: " + str(cmd),0)
            if fatal==1:raise OSError
        elif (error > 0):    
            my_print("MY_SYSTEM: System call is mislukt met errorcode: " + str(error),0)
            my_print("MY_SYSTEM: System call: " + str(cmd),0)
            if (fatal==1):raise OSError
    except:
        my_print("MY_SYSTEM: os.system is mislukt " + str(sys.exc_info()[0]),fatal)


def my_removefile(filename,my_sys_ldebug=0,fatal=0):
    """
    Verwijdert een file.
    Argumenten: pad naar file en evt een fatal-optie
    """
    #Controleer of opgegeven pad een bestaande file is.
    if os.path.isfile(filename):
        try:
            if (my_sys_ldebug):
                print("MY_REMOVEFILE: Deleting file: "+filename)
            os.remove(filename)
        except:
            my_print("MY_REMOVEFILE: Bestand " + str(filename) + " kan niet verwijderd worden.",fatal)
    else:
        my_print("MY_REMOVEFILE: Bestand " + str(filename) + " bestaat niet.",fatal)


   
def my_readfile(filename,my_sys_ldebug=0,fatal=1):
    """
    Leest een file in mbv readlines en geeft als return een list.
    Elke regel in de file is een string in de list.
    Argumenten: filename en evt een fatal-optie
    """
    # Controleer of opgegeven pad een bestaande file is.
    if os.path.isfile(filename):
        try:
            fpin = open(filename,"rU")
            if (my_sys_ldebug==1):
                print("MY_READFILE: Bestand " + filename + " is geopend.")
            list = fpin.readlines()
            fpin.close()
            del fpin
            if (my_sys_ldebug==1):
                print("MY_READFILE: Bestand " + filename\
                      +" is succesvol ingelezen met " + str(len(list))+" regels.")
            if (len(list)==0):
                my_print("MY_READFILE: U heeft een lege file ingelezen: " + str(filename),0)
            for regel in range(0,len(list)):
                # Als string eindigt met \n wordt deze verwijderd.
                if (list[regel][len(list[regel])-1:len(list[regel])]=="\n"):
                    list[regel]= list[regel][0:len(list[regel])-1]
            return list
        except:
            my_print("MY_READFILE: File " + str(filename) + " can not be read.",fatal)
    else:
        my_print("MY_READFILE: File " + str(filename) + " does not exist.",fatal)

def my_readfile2(filename,my_sys_ldebug=0,fatal=1):
    """
    Reads a file using read and returns a string
    Arguments: filename and optional fatal and my_sys_ldebug
    """
    # Check whether path is an existing file
    if os.path.isfile(filename):
        try:
            fpin = open(filename,"r")
            string = fpin.read()
            fpin.close()
            del fpin
            if (my_sys_ldebug==1):
                print("MY_READFILE: File ",filename," has ",len(string)," characters.")
            return string
        except:
            my_print("MY_READFILE2: File " + str(filename) + " can not be read.",fatal)
    else:
        my_print("MY_READFILE2: File " + str(filename) + " does not exist.",fatal)


def my_getcolumn(filename,collist,sep=None,my_sys_ldebug=0,fatal=1):
    """
    Selecteert uit een bestand alleen de gegevens van opgegeven kolommen 
    en geeft als return een list. Elementen van deze list zijn op hun beurt weer lists
    met daarin voor elke regel uit het bestand de gegevens uit de opgegeven kolommen.
    Argumenten: filename, list met kolomnummers, seperator (default = whitespace) 
    en evt een debug- en fatal-optie  NB: kolomnummers beginnen bij 0!
    """
    # Controleer eerst of opgegeven lijst met kolomnummers gevuld is.
    try:
        if len(collist)==0:
            my_print("MY_GETCOLUMN: U heeft geen kolomnummers opgegeven.",1)

        # Controleer daarna of opgegeven kolomnummers getallen zijn.        
        # Converteer colist van strings naar integers.
        try:
            icollist = []
            for col in range(0,len(collist)):
                icollist.append(int(collist[col]))
            if (my_sys_ldebug==1):
                print("MY_GETCOLUMN: Conversie van kolomnummers uit "\
                      +str(collist)+" naar integers is succesvol uitgevoerd.")
        except:
            my_print("MY_GETCOLUMN: Conversie van kolomnummers uit "\
                  + str(collist) + " naar integers is mislukt",fatal)

        # Bestand wordt ingelezen met my_readfile.
        if (my_sys_ldebug==1):
            print("MY_GETCOLUMN: Bestand "+filename+" wordt ingelezen.")
        list=my_readfile(filename,my_sys_ldebug,fatal)
        newlist=[]
        # Zoek voor elke regel in het bestand
        # de elementen uit de kolommen met de opgegeven kolomnummers. 
        for item in range(0,len(list)):
            fields = list[item].split(sep)
            hulplist=[]    
            for col in range(0,len(icollist)):
                # Controleer of opgegeven kolomnummers voorkomen in het bestand
                if ((icollist[col] > -1) and (icollist[col] < len(fields))):
                    # Plaats de gevonden elementen in een hulplist
                    hulplist.append(fields[icollist[col]])
                else:
                    my_print("MY_GETCOLUMN: Kolomnummer " + str(icollist[col]) + " is niet correct.",0)
                    my_print("MY_GETCOLUMN: In regel "+str(item+1)+" in bestand "+filename+\
                          " zijn "+str(len(fields))+" kolommen gevonden.",fatal)

            # Voeg de hulplisttoe aan de uiteindelijke list
            # en ga verder met de volgende regel in het bestand.
            newlist.append(hulplist)
        return newlist
    except:
        my_print("MY_GETCOLUMN: Functie my_getcolumn is mislukt met errorcode: " +\
                 str(sys.exc_info()[0]),fatal)

          
def my_getcolumnnumbers(filename,columnnames, sep=None,my_sys_ldebug=0,fatal=1):
    """
    Vergelijkt kolomnamen(case-insensitive) met text uit de header van een file
    en geef de bijpassende kolomnummers terug als list.
    Argumenten: filename, list met kolomnamen, seperator (default = whitespace) 
    en evt een debug- en fatal-optie. NB: kolomnummers beginnen bij 0!
    """
    # Controleer of opgegeven pad een bestaande file is
    if os.path.isfile(filename):
        try:
            numberlist = []
            fpin = open(filename,"r")
            if (my_sys_ldebug==1):
                print("MY_GETCOLUMNNUMBERS: Bestand "+filename+" is geopend.")
            # Lees alleen de eerste regel met de header.
            list = fpin.readline()
            # Remove \n van de ingelezen regel.
            list = list[0:len(list)-1]
            fpin.close()
            del fpin
            hulplist = list.split(sep)
            if (my_sys_ldebug==1):
                print("MY_GETCOLUMNNUMBERS: In de header van bestand  "\
                      +filename+" zijn "+str(len(hulplist))+" kolommen geteld.")
            #Vergelijk opgegeven kolomnamen een voor een met elementen uit de headerlist.
            for name in range(0,len(columnnames)):
                for item in range(0,len(hulplist)):
                    if (hulplist[item].upper()==columnnames[name].upper()):
                        numberlist.append(item)
                        #Als kolomnummer gevonden is, stap uit for loop
                        break
                #Als for-loop toch doorloopt tot einde, is kolomnummer niet gevonden.    
                else:
                    my_print("MY_GETCOLUMNNUMBERS: Kolom " +str(columnnames[name])\
                          + " komt niet voor in de header van "+str(filename),0)
        #Geef opgegeven kolomnamen weer als kolomnummers in een list.
            if (my_sys_ldebug==1):
                print("Laatst vergeleken kolomnaam is " +columnnames[name])
            return numberlist           
        except:
            my_print("MY_GETCOLUMNNUMBERS: Functie my_getcolumnumbers is mislukt\
                      met errorcode: " + str(sys.exc_info()[0]),fatal)
    else:
        my_print("MY_GETCOLUMNNUMBERS: Bestand " ,filename," bestaat niet.",fatal)

def my_get_environ_var(name,fatal=1):
    try:
        result = os.environ.get(name)
        if result is None:
            my_print("\nEnvironment variable " + str(name) + " must be set before starting the programme.",fatal)
        else:
            return result
    except:
        my_print("MY_GET_ENVIRON_VAR: Functie my_get_environ_var is mislukt met errorcode: " +\
                 str(sys.exc_info()[0]),fatal)


#Source code om de functies in deze module te testen.
           
if (__name__ == "__main__"):
    test_my_rmdir = 0
    test_my_copyfile = 1
    test_my_copytree = 0
    test_my_readfile = 0
    ldebug = 1

    # Test voor de my_rmdir
    # Vul hier een eigen directory in, die als test kan worden (aangemaakt/gekopieerd en) verwijderd
    if (test_my_rmdir): 
        if not os.path.isdir("d:\\temp\\applics"):
            os.mkdir("d:\\temp\\applics")
            os.system("xcopy /A /S d:\\data\\applics\\*.* d:\\temp\\applics")
        my_rmdir("d:\\temp\\applics")

    # Test voor de my_copyfile
    if (test_my_copyfile):
        file = "d:\\temp\\qq.txt"
        if not os.path.isfile(file):
            fp = open(file,"w")
            fp.write("lala")
            fp.close()
            del fp
        dir = os.path.split(file)[0]+"\\martine"
        if not os.path.isdir(dir):
            os.mkdir(dir) 
        my_copyfile(file,file+"new",1)
        # Nu checken op een destinatie als directory
        my_copyfile(file,dir,1)

    # Test voor de my_copyfile
    if (test_my_copytree):
        # Nu checken op een destinatie als directory
        my_rmdir("d:\\temp\\martine")
        my_copydir("d:\\temp\\martien","d:\\temp\\martine",1)
        
    # Test voor my_readfile
    if (test_my_readfile):
        # Test of teken ; aan het eind van laatste regel niet wordt verwijderd
        adres = "D:\\Programs\\GISMO\\datadirs\\countries\\econ\\"
        A = my_readfile(adres + "EconomicBase.dat")
        for line in range (220,len(A)):
            fields = A[line].split()
            for value in fields:
                print("Value: ",value)

