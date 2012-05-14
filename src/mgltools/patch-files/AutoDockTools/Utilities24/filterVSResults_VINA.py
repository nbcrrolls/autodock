#!/usr/bin/env python



#############################################################################
#
# Author: Stefano FORLI
#
# Copyright: Stefano Forli, TSRI 2011
#
# v.0.2
#############################################################################

from sys import argv
import getopt
import re
from glob import glob
import os
from AutoDockTools.HelperFunctionsN3P import pathToList, getLines, percent


# TODO
# -test "-d' mode
# -test "-F" mode
# -test -o mode
# -IDEA: support option for a config file? (list of interaction wewant and we dont"
# TODO remove the references to cluster and so on 

# TODO
# add an automated log file creation 
# filter_vina -F filtered_A-THR156_A-THR159_A-MET162.log  -H C:LYS68 -O -->
#                   filtered_A-THR156_A-THR159_A-MET162_C-LYS68.log
# TODO
# TODO
# Test the wildcard! "any" PHE!
# CORRECT THE USAGE FUNCTION FOR VINA
# TODO
#
# MIGRATE TO THE GNU-COMPLIANT DOUBLE DASH
# http://www.doughellmann.com/PyMOTW/getopt/
#

# defaults
mode = "1"   # default pose
ebest = -999.    # energy
eworst = -3.  
cworst = 1.     # cluster poses
cbest = 999.
pworst = 1.     # cluster %
pbest = 100.
lworst = 0.     # ligand efficiency
lbest = -99
DEBUG = False
pattern = "*_VS.pdbqt" # default pattern for searching result files
do_filter = False
ecursive = False

def checkVSresult(lines):
    "parses and validates a PDBQT+ file (lines)"
    if lines[0].startswith("USER    ADVS_Vina_result>"):
        return True
    else:
        return False

def getResultCount(lines):
    return int(lines[4].split("ADVina_results>")[1].strip())

def getHistogram(lines):
    return lines[5].split("ADVina_histogram>")[1].strip().split(",")

def setKwMode(mode = "1"):
    # USER    ADVina_pose1> -10.200,	-0.309
    # mode  = "1", "2", "..."
    if mode == "any":
        kw = "ADVina_pose."
    else: 
        kw = "ADVina_pose"+mode
    return kw

def getRawEnergyData(lines, mode = "1"):
    # mode  = "1", "2", "...", "any" (POSES)
    # format:  { "e" : [ float(e)] , "leff" : [float(leff)] }
    kw = setKwMode(mode)+">"
    result =  { "e"     : [], 
                "leff"  : [] }
    for l in lines:
        if re.search(kw, l):
            l = l.split(">", 1)[1]
            e, leff = l.split(",")    
            result['e'].append(float(e))
            result["leff"].append(float(leff))
            if not mode == "any":
                break
    return result


def getLigInteractions(lines, mode = "1"):
    kw = setKwMode(mode) # 
    kw += "_" # TODO test this!
    interactions =  {"vdw" : [],  # these keys must match the 
                    "hba" : [],   # tags used in writing the PDBQT+
                    "hbd" : [],
                    "ppi" : [],
                    "tpi" : [],
                    "mtl" : [],
                    }
    for l in lines:
        if re.search(kw, l):
            for itype in interactions.keys():
                if (itype == "ppi") or (itype == "tpi"):
                    sep = ";"
                else:
                    sep = ","
                if itype in l:
                    l = l.split(itype+">", 1)[1]
                    l = l.split(sep)
                    for i in l:
                        interactions[itype].append(i.strip())
    if DEBUG:
        print "=== Extracted interactions ==="
        for i in interactions:
            print i, interactions[i]

    return interactions

def filterLigandProperties(ligdata, eworst, ebest, lworst, lbest):
    """
    INPUT : ligdata format:  { "e" : [ float(e)] , "leff" : [float(leff)] }
    OUTPUT: Bool(ligdata satisfy requested properties)
    NOTES : every property is a list to handle requests where any ligand conformation can
            satisfy the properties.
    """
    KEEP = False
    for x in ligdata['e']:
        if x <= eworst and x >= ebest:
            KEEP = True
            break
    if not KEEP:
        if DEBUG: print "Energy violated", ligdata['e']
        return False
    else:   
        KEEP = False
    for x in ligdata['leff']:
        if x <= lworst and x >= lbest:
            KEEP = True
            break
    if not KEEP:
        if DEBUG: print "Ligand efficiency violated", ligdata['leff']
        return False
    else:   
        KEEP = True
    # trunkated here
    return True


def matchThisInteraction(entry1, entry2):
    "checks if entry1 is equal or a subset of entry2"
    # entry1 = ":TYR12:"
    # entry2 = "A:TYR12:O"
    #print entry1, entry2
    parts = entry1.count(":")
    if entry1.count(":") == 2:      # vdw, hb, metal  (D:DA17:N1)
        chain1, res1, at1 = entry1.split(":")
    elif entry1.count(":") == 1:    # pi interaction  (A:MG396:MG)
        chain1, res1 = entry1.split(":")
        at1 = ""
    if entry2.count(":") == 2:      # vdw, hb, metal  (D:DA17:N1)
        chain2, res2, at2 = entry2.split(":")
    elif entry2.count(":") == 1:    # pi interaction  (A:MG396:MG)
        chain2, res2 = entry2.split(":")
        at2 = ""

    if DEBUG: 
        pass
        #print "chain, res, at: comparing..."
        #print "|%s:%s:%s|" %(chain1, res1, at1),
        #print "|%s:%s:%s|" %(chain2, res2, at2),

    if not chain1 or chain1 == chain2:
        if not res1 or res1 == res2:
            if not at1 or at1 == at2:
                if DEBUG: print "found!"
                return True
    return False




# interaction filters functions and sub-functions
def filterVdw(liginteractions, requested):
    found = 0
    for c in requested: # c = "[-]chain:resNameNum:atom"
        if c[0] == "-":
            c = c[1:]
            for atom in liginteractions:
                if matchThisInteraction(c, atom): 
                    if DEBUG: print "FAIL [vdw]: found unwanted interaction", c
                    return False 
            found += 1
        else:
            for atom in liginteractions:
                if matchThisInteraction(c, atom): 
                    found += 1
                    break
    if not found == len(requested):
        if DEBUG: print "FAIL [vdw]", found, len(requested)
        return False
    else: return True
    
def filterHba(liginteractions, requested):
    found = 0
    for c in requested:
        if c[0] == "-":
            c = c[1:]
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[1]): 
                    if DEBUG: print "FAIL [hba]: found unwanted interaction", c
                    return False 
            found += 1
        else:
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[1]):
                    found += 1
                    break
    if not found == len(requested):
        if DEBUG: print "failed hba", found, len(requested)
        return False
    else: return True

def filterHbd(liginteractions, requested):
    found = 0
    for c in requested:
        if c[0] == "-":
            c = c[1:]
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[1]):
                    if DEBUG: print "FAIL [hbd]: found unwanted interaction", c
                    return False 
            found += 1
        else:
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[1]):
                    found += 1
                    break
    if not found == len(requested):
        if DEBUG: print "failed hbd"
        return False
    else: return True

def filterHb(liginteractions, requested):
    found = 0
    for c in requested:
        if c[0] == "-":
            c = c[1:]
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[1]):
                    if DEBUG: print "FAIL [hb]: found unwanted interaction", c
                    return False 
            found += 1
        else:
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[1]):
                    found += 1
                    break
    if not found == len(requested):
        if DEBUG: print "failed hb (generic)"
        return False
    else: return True

def filterPpi(liginteractions, requested):
    found = 0
    for c in requested:
        if c[0] == "-":
            c = c[1:]
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[0]):
                    if DEBUG: print "FAIL [ppi]: found unwanted interaction", c
                    return False 
            found += 1
        else:
            for atom in liginteractions:
                if (matchThisInteraction(c,atom.split("~~")[0]) == desired):
                    found += 1
                    break
    if not found == len(requested):
        if DEBUG: print "failed ppi", found, len(requested)
        return False
    else: return True

def filterTpi(liginteractions, requested):
    found = 0
    for c in requested:
        if c[0] == "-":
            desired = -1
            c = c[1:]
        else:
            desired = 1
        for atom in liginteractions:
            if (matchThisInteraction(c,atom.split("~~")[0]) == desired):
                found += 1
    if not found == len(requested):
        if DEBUG: print "failed tpi", found, len(requested)
        return False
    else: return True

def filterPi(liginteractions, requested):
    found = 0
    for c in requested:
        if c[0] == "-":
            c = c[1:]
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[0]):
                    if DEBUG: print "FAIL [pi]: found unwanted interaction", c
                    return False 
            found += 1
        else:
            for atom in liginteractions:
                if (matchThisInteraction(c,atom.split("~~")[0]) == desired):
                    found += 1
                    break
    if not found == len(requested):
        if DEBUG: print "!!!FAIL [pi] ", found, len(requested)
        return False
    else: return True

def filterMtl(liginteractions, requested):
    found = 0
    for c in requested:
        if c[0] == "-":
            c = c[1:]
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[1]): 
                    if DEBUG: print "FAIL [mtl]: found unwanted interaction", c
                    return False 
            found += 1
        else:
            for atom in liginteractions:
                if matchThisInteraction(c, atom.split("~~")[1]):
                    found += 1
                    break
    if not found == len(requested):
        if DEBUG: print "failed metal", found, len(requested)
        return False
    else: return True


def countHb(liginteractions, requested):
    if requested == 0: return True
    if DEBUG: print "HB COUNT: %d | %d [ requested, found ]" %( requested, len(liginteractions)),
    if requested <= len(liginteractions):
        if DEBUG: print "=> OK"
        return True
    if DEBUG: print "=> Fail"
    return False

def filterLigandInteractions(liginteractions, int_filter):
    # van der Waals 
    if not filterVdw(liginteractions['vdw'], int_filter['vdw']):
        return False
    # hydrogen bond acceptor 
    if not filterHba(liginteractions['hba'], int_filter['hba']):
        return False
    # hydrogen bond donor 
    if not filterHbd(liginteractions['hbd'], int_filter['hbd']):
        return False
    # hydrogen bond (any kind)
    if not filterHb(liginteractions['hba']+liginteractions['hbd'], int_filter['hb']):
        return False
    # pi interaction (pi-pi stacking)
    if not filterPpi(liginteractions['ppi'], int_filter['ppi']):
        return False
    # pi interaction (T-stacking)
    if not filterPpi(liginteractions['tpi'], int_filter['tpi']):
        return False
    # pi interaction (any kind)
    if not filterPi(liginteractions['tpi']+liginteractions['ppi'], int_filter['pi']):
        return False
    # Metal
    if not filterMtl(liginteractions['mtl'], int_filter['mtl']):
        return False
    # number of hb donors (ligand side)
    if not countHb(liginteractions['hbd'], int_filter['hbd_count']):
        return False
    # number of hb acceptors (ligand side)
    if not countHb(liginteractions['hba'], int_filter['hba_count']):
        return False
    # number of hb (ligand side)
    if not countHb(liginteractions['hbd']+liginteractions['hba'], int_filter['hb_count']):
        return False
    # if we got here...
    if DEBUG: print "ligand passed the test! ==> ",
    return True


######################################################################################

def usage():
    print """
USAGE
   %s  -f filename | -d directory | -F   [ OPTIONS ]
   -f filename     ligand PDBQT file to filter
   -d directory    directory containing PDBQT files to filter
   -F filename     file containing the list of files to be filtered
                   (i.e. generated with '-o' option)

OPTIONS
   -o filename     specify the output list file where filtered results saved
                   (i.e. to be used later with '-F' option
   -m mode         specify which pose must be considered the result [default: 1]
   -R              (with '-d') recursive search for ligand files in sub-directories
   -s pattern      (with '-d') search for ligand filenames containing 'pattern' [default '_VS.pdbqt']
PROPERTY FILTERS
   -e best:worst   energy cut-off [ default -999.0::-3.0 Kcal/mol ]
   -l best:worst   ligand efficiency filter [ default -99:0 Kcal ]
   -c best:worst   poses count in result cluster [ default 999:1 poses ]
   -p best:worst   percent poses in result cluster [ default 100:1 %% ]

    Filters can be specified as single values or as ranges:
          -e -7             --> energy of -7.0 Kcal/mol or better
          -l -0.6:-0.2      --> results with ligand efficiency between -0.6 and -0.2
          -c 5:100          --> cluster population between 5 and 100 poses
          -e :-5            --> energy no better than -5.0 Kcal/mol

INTERACTION FILTERS
    -V chain:res:atm    van der Waals interactions
    -A chain:res:atm    hydrogen bond (ligand as acceptor)
    -D chain:res:atm    hydrogen bond (ligand as donor)
    -H chain:res:atm    hydrogen bond (ligand as either donor or acceptor)
    -M chain:res:atm    metal coordination                
    -P chain:res        pi-interactions (parallel stacking) 
    -T chain:res        pi-interactions (T-stacking)     
    -L chain:res        any lone pair pi interaction
    -N number           minimum number of hb donated  (ligand as donor)
    -X number           minimum number of hb accepted (ligand as acceptor)
    -S number           minimum number of hb formed   (any type)

    Interaction keys can be specified as chain:residue, or separately, es.:
          -V B:LEU63:CG      --> contact interaction with carbon CG of Leu63 on chain B
          -A C:              --> hb acceptor interaction with any residue on chain C
          -D A:ASP25,B:ASP25 --> hb donor with both ASP25 on chains A and B.
          -H :LYS            --> any hydrogen bond interaction with any LYS
          -V -B:THR276       --> avoid contact with with THR276

       NOTES: 
         -residue and chain names are *cAsE sEnSiTiVe*
         -unwanted interactions can be specified by adding a '-' before the residue.

EXAMPLE
    %s -d . -m any -e :-10 -l -0.35:-0.6 -p 51 -V B:LEU63

    Filter all ligands present in the current directory for which either LC or LE have :
          -energy of -10.0 Kcal/mol or worst
          -ligand efficiency between -0.35 and -0.6 Kcal/mol
          -cluster population of 51%% minimum
          -close contact with Leu63 on chain B
""" % (argv[0], argv[0])


try:
    options, var = getopt.getopt(argv[1:], 'f:F:d:m:e:l:L:c:p:vV:A:D:H:P:T:M:o:Rs:N:X:S:')
except getopt.GetoptError, err:
    usage()
    print str(err), "\n"
    exit(2)

# user input validation
opts = {} # 
for o, a in options:
    # TODO TODO TODO  IT WAS VERIFIED TO BE *EXTREMELY* HELPFUL SO IT HAS TO BE DONE
    # Experimental: testing multiple instances
    # of the same option! (for interactions)
    #print "processing",o,a
    if o in opts.keys():
        opts[o] += ","+a
    else:
        opts[o] = a
    # TODO TODO TODO test this
    #opts[o] = a
if "-v" in opts:
    DEBUG = True
if (not "-f" in opts) and (not "-d" in opts) and (not '-F' in opts) or ('-f' in opts and '-d' in opts):
    print "\n\n ERROR: Either file (-f) or directory (-d) input must be specified! ###\n"
    usage()
    exit(1)
try:
    if '-m' in opts: 
        mode = opts['-m']
    if '-e' in opts:
        try:
            value = opts['-e']
            if ":" in value:
                if value[0] == ":" :
                    ebest = float(value[1:])
                else:
                    tmp = value.split(":")
                    eworst = float(tmp[0])
                    ebest = float(tmp[1])
                    if eworst < ebest:
                        print "WARNING: Energy min. (%2.2f) is better than max. (%2.2f)!" % (eworst, ebest)
            else:
                eworst = float(value)
        except:
            print "ERROR: problem parsing energy filter values '%s'" % opts['-e']
            usage()
            exit(1)
    if '-l' in opts:
        try:
            value = opts['-l']
            if ":" in value:
                if value[0] == ":" :
                    lbest = float(value[1:])
                else:
                    tmp = value.split(":")
                    lworst = float(tmp[0])
                    lbest = float(tmp[1])
                    if lworst > lbest:
                        print "WARNING: Ligand efficiency min. (%2.2f) is better than max. (%2.2f) " % (lworst, lbest)
            else:
                lworst = float(value)
        except:
            print "ERROR: problem parsing lig.efficiency filter values '%s'" % opts['-l']
            usage()
            exit(1)

    if '-c' in opts:
        try:
            value = opts['-c']
            if ":" in value:
                if value[0] == ":" :
                    cbest = float(value[1:])
                else:
                    tmp = value.split(":")
                    cworst = float(tmp[0])
                    cbest = float(tmp[1])
                    if cworst > cbest:
                        print "WARNING: Cluster size min. (%2.2f) is bigger than max. (%2.2f) " % (cworst, cbest)
            else:
                cworst = float(value)
        except:
            print "ERROR: problem parsing cluster filter values '%s'" % opts['-c']
            usage()
            exit(1)
        
    if '-p' in opts:
        try:
            value = opts['-p']
            if ":" in value:
                if value[0] == ":" :
                    pbest = float(value[1:])
                else:
                    tmp = value.split(":")
                    pworst = float(tmp[0])
                    pbest = float(tmp[1])
                    if cworst > cbest:
                        print "WARNING: Cluster size %% min. (%2.2f) is bigger than max. (%2.2f) " % (pworst, pbest)
            else:
                pworst = float(value)
        except:
            print "ERROR: problem parsing cluster%% filter values '%s'" % opts['-p']
            usage()
            exit(1)

    if '-s' in opts:
        pattern = opts['-s'] 
    if DEBUG: print "Pattern : |%s|" % pattern


    if '-f' in opts:
        files = [ opts['-f'] ] 
    elif '-d' in opts:
        if '-R' in opts:
            recursive = True
        files = pathToList( path = opts['-d'], recursive = recursive, pattern = pattern )
        if not files:
            print "WARNING: no VS results in the dir:", opts['-d'] 
            exit(1)
    elif '-F' in opts:
        try:
            files = getLines(opts['-F'], doStrip = True)
            print "%d files from list file %s" % (len(files), opts['-F'])
        except:
            print "ERROR opening the input list: " , opts['-F']
            exit(1)
    else:
        print "ERROR: Either file (-f), directory (-d) or list file (-F) input options must be specified"
        usage()
        exit(1)

    if '-R' in opts and ('-f' in opts or '-F' in opts):
        print "WARNING: recursive option ('-R') ignored in file/list mode."

    # interaction filters are initialized
    int_filter = { 'vdw' : [],
                   'hba': [],
                   'hbd': [],
                   'hb': [],
                   'pi'  : [],
                   'ppi' : [],
                   'tpi' : [],
                   'mtl':[],
                   'hbd_count' : 0,
                   'hba_count' : 0,
                   'hb_count' : 0 }
    if '-V' in opts:
        flt_contacts = opts['-V'].split(",")
        int_filter['vdw'] += filter(lambda x: x, flt_contacts) # purge empty values
        do_filter = True
    if '-A' in opts:
        flt_hacc = opts['-A'].split(",")
        int_filter['hba'] += filter(lambda x: x, flt_hacc)
        do_filter = True
    if '-D' in opts:
        flt_hdon = opts['-D'].split(",")
        int_filter['hbd'] += filter(lambda x: x, flt_hdon)
        do_filter = True
    if '-H' in opts:
        flt_hany = opts['-H'].split(",")
        int_filter['hb'] += filter(lambda x: x, flt_hany)
        do_filter = True
    if '-P' in opts:
        flt_ppi = opts['-P'].split(",")
        int_filter['ppi'] += filter(lambda x: x, flt_ppi)
        do_filter = True
    if '-T' in opts:
        flt_tpi = opts['-T'].split(",")
        int_filter['tpi'] += filter(lambda x: x, flt_tpi)
        do_filter = True
    if '-L' in opts:
        flt_pi = opts['-L'].split(",")
        int_filter['pi'] += filter(lambda x: x, flt_pi)
        do_filter = True
    if '-M' in opts:
        flt_metal = opts['-M'].split(",")
        int_filter['mtl'] += filter(lambda x: x, flt_metal)
        do_filter = True
    if '-N' in opts:
        int_filter['hbd_count'] = int(opts['-N'])
        do_filter = True
    if '-X' in opts:
        int_filter['hba_count'] = int(opts['-X'])
        do_filter = True
    if '-S' in opts:
        int_filter['hb_count'] = int(opts['-S'])
        do_filter = True
    if DEBUG:
        print "############## INTERACTION FILTERS ######################"
        for i in int_filter.keys():
            print "# ", i, int_filter[i]
        print "#######################################################\n"
except getopt.GetoptError, err:
    usage()
    print str(err), "\n"
    exit(2)

# /end user input validation
if DEBUG:
    print "############## PROPERTIES FILTER ######################"
    print "#                  min     MAX"
    print "#  --------------------------------"
    print "#  Energy      :", eworst,",", ebest
    print "#  Ligand eff. :", lworst,",", lbest
    print "#  Cluster #   :", cworst,",", cbest
    print "#  Cluster %   :", pworst,",", pbest
    print "#  Mode        :", mode 
    print "#######################################################\n"

if '-o' in opts:
    try: 
        logfile = opts['-o']
        if os.path.isfile(logfile):
            renamed = logfile+"_OLD"
            print "WARNING: log file '%s' exists. ( renamed => '%s' )" % (logfile, renamed)
            os.rename(logfile, renamed)
        output = open(logfile, 'w')
    except: 
        print "ERROR opening the output list file :" % opts['-o']
        exit(1)
    

tot = len(files)

c = 0
if DEBUG: print "\nSTARTING FILTERING...=================="
for f in files:
    #print "processing:", f
    try:
    #if True:
        ligand = getLines(f)
        if checkVSresult(ligand):
            if getResultCount(ligand) == 1 and not mode == "1": # TODO fix this!
                if DEBUG: print "[ single result => POSE 1 ]",
                current_mode = "1"
            else:
                current_mode = mode
            ligdata = getRawEnergyData(ligand, mode = current_mode)
            passed = False
            if filterLigandProperties(ligdata, eworst, ebest, lworst, lbest):
                if do_filter:
                    liginteractions = getLigInteractions(ligand, mode = current_mode)
                    if filterLigandInteractions(liginteractions, int_filter):
                        passed = True
                        c+=1
                else:
                    passed = True
                    c+=1
            if passed:
                if '-o' in opts:
                    try:
                        output.write(f.strip()+"\n")
                    except:
                        print "ERROR writing the output file list"
                        exit(2)
                else:
                    print f
        else:
            print "ligand is not a good VS", f
    except:
    #else:
        print "# ERROR with file : %s" % f
        exit(1)
print "%d filtered => %d passed" % (len(files),c)


