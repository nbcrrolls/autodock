#!/usr/bin/env python



#############################################################################
#
# Author: Stefano FORLI
#
# Copyright: Stefano Forli, TSRI 2011
#
# v.0.4
#############################################################################

from sys import argv
import getopt
import re
from glob import glob
import os
from AutoDockTools.HelperFunctionsN3P import pathToList, getLines, percent
#TODO from AutoDockTools.HelperFunctionsN3P import pathToList, getLines, percent, matchThisInteraction
# TODO 

version = '0.4beta'

# TODO 
# -IDEA: support option for a config file? (list of interaction wewant and we dont"
# TODO migrate 
# TODO

# defaults
mode = "lc"   # default pose
ebest = -999.    # energy
eworst = -3.  
cworst = 1.     # cluster poses
cbest = 999.
pworst = 1.     # cluster %
pbest = 100.
lworst = 0.     # ligand efficiency
lbest = -99
DEBUG = False
pattern = "_VS.pdbqt" # default pattern for searching result files
do_filter = False
recursive = False

def checkVSresult(lines):
    "parses and validates a PDBQT+ file (lines)"
    if lines[0].startswith("USER    ADVS_result>"):
        return True
    else:
        return False

def getRecName(lines):
    return lines[1].split("AD_rec>")[1].strip()

def getRawRunsData(lines): # TODO USELESS?
    return lines[2].split("AD_runs,rmstol,tot_clusters>")[1].strip().split(",")

def getDlgList(lines):
    return lines[3].split("AD_dlg_list>")[1].strip().split(",")

def getResultCount(lines):
    return int(lines[4].split("AD_results>")[1].strip())

def getHistogram(lines):
    return lines[5].split("AD_histogram>")[1].strip().split(",")

def setKwMode(mode = "lc"):
    # mode  = "le", "lc", "any"
    if mode == "le":
        kw = "AD_LE"
    elif mode == "lc":
        kw = "AD_LC"
    elif mode == "any":
        kw = "AD_L."
    return kw

def getRawEnergyData(lines, mode = "lc"):
    # mode  = "le", "lc", "any"
    # format:  { "e" : [ float(e)] , "leff" : [float(leff)], "c_size" : [int(c_size)], "c_pc" : [float(c_pc)] }
    kw = setKwMode(mode)+">"
    result =  { "e"     : [], 
                "leff"  : [],
                "c_size": [],
                "c_pc"  : [] }
    for l in lines:
        if re.search(kw, l):
            l = l.split(">", 1)[1]
            e, leff, c_size, c_pc = l.split(",")    
            result['e'].append(float(e))
            result["leff"].append(float(leff))
            result["c_size"].append(int(c_size))
            result["c_pc"].append(float(c_pc))
            if not mode == "any":
                break
    return result


def getLigInteractions(lines, mode = "lc"):
    kw = setKwMode(mode) # 
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

def filterLigandProperties(ligdata, eworst, ebest, lworst, lbest, cworst,cbest, pworst,pbest):
    """
    INPUT : ligdata format:  { "e" : [ float(e)] , "leff" : [float(leff)], "c_size" : [int(c_size)], "c_pc" : [float(c_pc)] }
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
        KEEP = False
    for x in ligdata['c_size']:
        if x >= cworst and x <= cbest:
            KEEP = True
            break
    if not KEEP:
        if DEBUG: print "Cluster size violated",ligdata['c_size']
        return False
    else:   
        KEEP = False
    for x in ligdata['c_pc']:
        if x >= pworst  and x <= pbest:
            return True
    if DEBUG: print "Cluster percent violated", ligdata['c_pc']
    return False


def matchThisInteraction(entry1, entry2, DEBUG=False):
    "checks if entry1 is equal or a subset of entry2"
    # entry1 = ":TYR12:"
    # entry2 = "A:TYR12:O"
    if DEBUG:
        print "E1:", entry1
        print "E2:", entry2
    parts = entry1.count(":")
    if parts == 2:      # vdw, hb, metal  (D:DA17:N1)
        chain1, res1, at1 = entry1.split(":")
    elif parts == 1:    # pi interaction  (A:MG396:MG)
        chain1, res1 = entry1.split(":")
        at1 = ""

    chain2, res2, at2 = entry2.split(":")
    #chain2, res2 = entry2.split(":")

    if DEBUG: 
        print "chain, res, at: comparing..."
        print "|%s:%s:%s|" %(chain1, res1, at1),
        print "|%s:%s:%s|" %(chain2, res2, at2),

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
    print """                                    (AutoDockTools)

NAME
        AutoDock Results Filter v.%s
        
SYNOPSIS
        %s  -f filename | -d directory | -F  [OPTIONS] [FILTERS]

DESCRIPTION
        Filter AutoDock result files (PDBQT+) by properties (energy, 
        ligand efficiency...) and interactions calculated with the 
        receptor structure. 
        Input files can be defined by specified a single filename, 
        a directory to scan or a list file containing the paths to 
        files to scan.
 
        -f filename     single ligand PDBQT+ file to filter
        -d directory    directory containing PDBQT+ files to filter
        -F filename     file containing PDBQT+ filenames
 
OPTIONS
        -o filename     specify the output list file where filtered
                        results saved
        -m mode         specify which pose must be considered the result:
                         lc   (lowest energy in the largest cluster) [default]
                         le   (lowest energy pose is selected)
                         any  (either lowest energy or largest cluster)
        -R              recursive search for ligand files in sub-directories
                            [only with '-d']
        -s pattern      search for ligand filenames containing 'pattern' 
                            [only with '-d', default '_VS.pdbqt']
 
PROPERTY FILTERS 
        -e min[:max]   energy range [ default -0.3:-999.0 Kcal/mol ]
        -l min[:max]   ligand efficiency range [ default 0:-99 Kcal ]
        -c min[:max]   poses count in result cluster [ default 1:999 poses ]
        -p min[:max]   percent poses in result cluster [ default 1:100 %% ]
 
        Filters can be specified as single values or as ranges:
          -e -7             --> energy of -7.0 Kcal/mol or better
          -l -0.2:-0.6      --> ligand efficiency between -0.6 and -0.2
          -c 5:100          --> cluster population between 5 and 100 poses
          -e :-5            --> energy no better than -5.0 Kcal/mol
                                ( [+oO,-0.5] )
 
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
 
        Interaction keys can be specified (or negated) as chain:residue:atom, 
        or separately, es.:
             -V B:LEU63:CG      --> contact interaction with carbon CG of 
                                    Leu63 on chain B
             -A C:              --> hb acceptor interaction with any residue
                                    on chain C
             -D A:ASP25,B:ASP25 --> hb donor with both ASP25 on chains A and B.
             -H :LYS            --> any hydrogen bond interaction with any LYS
             -V -B:THR276       --> *avoid* contact with with THR276
        Residue and chain names are *cAsE sEnSiTiVe*.
        Unwanted interactions can be specified by adding a '-' before the 
        residue string (i.e. -B:THR:276).

EXAMPLES  
        %s -d . -m any -e -10 -l -0.35:-0.6 -p 51 -V B:LEU63
  
        Filter all ligands present in the current directory for which either
        LC or LE satisfy the following filters:
             -energy of -10.0 Kcal/mol or worst
             -ligand efficiency between -0.35 and -0.6 Kcal/mol
             -cluster population of 51%% minimum
             -close contact with Leu63 on chain B

AUTHOR
        Written by Stefano Forli. Contributions by Alex L. Perryman.

REPORTING BUGS
        Please report bugs to:
        MGL BugZilla            http://mgldev.scripps.edu/bugs/
        MGL forum               http://mgl.scripps.edu/forum
        AutoDock mailing list   http://autodock.scripps.edu/mailing_list

COPYRIGHT
        Copyright (C) 2011 Stefano Forli, Molecular Graphics Lab, 
                     The Scripps Research Institute.
        GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
""" % (version, os.path.basename(argv[0]), os.path.basename(argv[0]))

try:
    options, var = getopt.getopt(argv[1:], 'f:F:d:m:e:l:L:c:p:vV:A:D:H:P:T:M:o:Rs:N:X:S:')
except getopt.GetoptError, err:
    usage()
    print str(err), "\n"
    exit(2)

# user input validation
opts = {} # 
for o, a in options:
    # TODO TODO TODO 
    # Experimental: testing multiple instances
    # of the same option! (for interactions)
    #print "processing",o,a
    #if o in opts.keys():
    #    opts[o] += ","+a
    #else:
    #    opts[o] = a
    # TODO TODO TODO test this
    opts[o] = a

if len(opts) == 0:
    usage()
    exit(0)

if "-v" in opts:
    DEBUG = True
if not "-f" in opts and not "-d" in opts or ('-f' in opts and '-d' in opts):
    print "\n\n ### Error! Either file (-f) or directory (-d) input must be specified! ###\n"
    usage()
    exit(1)

try:
    if '-R' in opts:
        recursive = True
    if '-m' in opts: 
        mode = opts['-m']
    if '-e' in opts:
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
    if '-l' in opts:
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
    if '-c' in opts:
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
    if '-p' in opts:
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
    
    if '-s' in opts:
        pattern =  opts['-s']
        if DEBUG: print "filename pattern search>", pattern
    if '-f' in opts:
        files = [ opts['-f'] ] 
    elif '-d' in opts:
        files = pathToList( path = opts['-d'], recursive = recursive, pattern = pattern )
        if not files:
            print "No VS results in the dir:", opts['-d'] 
            exit(1)
    elif '-F' in opts:
        try:
            files = getLines(opts['-L'])
        except:
            print "ERROR opening the input list: " % opts['-L']
            exit(1)
    else:
        print "ERROR: Either file (-f), directory (-d) or list file (-F) input options must be specified"
        usage()
        exit(1)
        
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
        print "############## REQUESTED INTERACTION FILTERS ######################"
        for i in int_filter.keys():
            print "# ", i, int_filter[i]
        print "###################################################################\n"
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
        output = open(opts['-o'], 'w')
    except:
        print "Error opening the output list file :" % opts['-o']
        exit(1)


problematic = None

for f in files:
    try:
    #if True:
        ligand = getLines(f)
        if checkVSresult(ligand):
            if getResultCount(ligand) == 1 and mode == "lc":
                if DEBUG: print "[ single result => LE ]",
                current_mode = "le"
            else:
                current_mode = mode
            ligdata = getRawEnergyData(ligand, mode = current_mode)
            passed = False
            if filterLigandProperties(ligdata, eworst, ebest, lworst, lbest, cworst, cbest, pworst, pbest):
                if DEBUG: print "passed properties filter |",
                if do_filter:
                    liginteractions = getLigInteractions(ligand, mode = current_mode)
                    # filter contacts:
                    if filterLigandInteractions(liginteractions, int_filter):
                        passed = True
                        #if DEBUG: print "passed interactions."
                else:
                    passed = True
            else:
                if DEBUG: print "Failed properties filter!"

            if passed:
                if '-o' in opts:
                    try:
                        output.write(f.strip()+"\n")
                    except:
                        print "ERROR writing the output file list [DISK FULL?]"
                        exit(2)
                else:
                    print f
        else:
            if DEBUG: print "ligand didn't pass", f
    except:
    #else:
        if DEBUG: print "# ERROR with file : %s" % f
        if not problematic:
            try:
                problematic = open("problematic.log", 'w')
            except:
                print "error in writing the 'problematic.log' file.\n Aborting..."
                exit(1)
        problematic.write(f+"\n")

if '-o' in opts:
    try:
        output.close()
    except:
        print "ERROR when closing log file '%s'." % (opts['-o'])
        exit(1)
if problematic:
    problematic.close()

exit(0)


