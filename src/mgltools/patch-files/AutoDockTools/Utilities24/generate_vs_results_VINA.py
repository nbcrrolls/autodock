#!/usr/bin/python

#############################################################################
#
# Author: Stefano FORLI
#
# Copyright: Stefano Forli, TSRI 2011
#
# v.0.2
#############################################################################




# TODO add the automatic scanning function for getting DLGs, identifying the ligands and so on...

########################################################################################
# specific helper functions
def doDebug(lig):
    if doInteractions:
        buffer = []
        try:
            fp = open(lig.ligName+"_DEB_contacts.pdb", 'w')
            fp.write("MODEL\nUSER   vdw contacts \n")
            for x in lig.results[0]['vdw_contacts']:
                fp.write(x+'\n')
            fp.write('ENDMDL \n')
            fp.close()

            fp = open(lig.ligName+"_DEB_hb_accept.pdb", 'w')
            fp.write("MODEL\nUSER   hb acceptors\n")
            for p in lig.results[0]['hb']['acceptors']:
                for x in p:
                    fp.write(x+'\n')
            fp.write('ENDMDL \n')
            fp = open(lig.ligName+"_DEB_hb_donor.pdb", 'w')
            fp.write("MODEL\nUSER   hb donors\n")
            for p in lig.results[0]['hb']['donors']:
                for x in p:
                    fp.write(x+'\n')
            fp.write('ENDMDL \n')
            fp.close()
        except:
            print "ERROR! writing the debug files"
            exit(2)

def usage():
    print """

USAGE    %s  [-f filename|-d directory] -r receptor [ OPTIONS ]


            -f filename     Vina Out PDBQT file(s) to process. Multiple files can be separated by commas:
                            "ligand1_out_pdbqt,ligand2_out.pdbqt,..."
            -d directory    directory containing Vina output files to process
            -r receptor     specify receptor filename (PDBQT)

OPTIONS
            -i              disable ligand-target interactions calculation 
            -l filename     write ligand results summary to tab separated filename
                            - if the file does not exist, it will be created
                            - if the file exists, results ligand results will be append
            -R              scan 'directory' recursively

            -m X            specify the maximum number of poses to get from the Vina results
                            (default: 1)
            -p pattern      process all the files containing the pattern [default : "*_out.pdbqt"]
            -v              set verbose mode On (DEBUG)

            -x X.XX         set the extra tolerance (Angstroms) for the hydrogen bond interaction

EXAMPLE
            %s -d . -r protein.pdbqt -t 1.0 

            Process all the Vina results in the current directory with 1.0A RMSD tolerance and use "protein.pdbqt"
            coordinates to calculate ligand-target interactions.
""" % (argv[0], argv[0])

#            -R              scan specified directory and all sub-directories for DLG(s)
#                            ( ignored with -f)

#===============================================
if __name__ == '__main__':
    from glob import glob  
    import os
    from sys import argv, exc_info
    import getopt
    from AutoDockTools.VsResultsGenerator import *
    from AutoDockTools.HelperFunctionsN3P import pathToList, writeList

    # defaults
    DEBUG = False
    doInteractions = True
    recursive = False
    rmsTol = 2.0
    mode = 1 # binding modes extracted from the result 
    pattern = "*_out.pdbqt"
    header = "#name\tenergy\tligand_efficiency\ttotal_poses\tfilename\n"
    suffix = "_Vina_VS"
    hbtol = 0.0

    try:
        options, var = getopt.getopt(argv[1:], 'f:d:r:o:ivm:l:p:Rx:')
        # -f  str       dlg single file
        # -d  str       directory
        # -r  str       receptor file
        # -R  (none)    recursive (find all Vina results in a given path matching patternt)
        # -o  str       output filename (override defult "ligname_VS.pdbqt")
        # -m  int       set the number of modes that need to be extracted (default 1)
        # -l  logfile   add the logfile info  TODO 
        # 
        # -x  float     extra tolerance for hbond.
        # -i            don't calculate interactions
        # -v  (none)    verbose (DEBUG MODE)

        # TODO very important!!!
        # TODO support the caching of the receptor? Why not?....
        # TODO ADD THE EXTENSION OPTION IN TEH LIST!!!

    except getopt.GetoptError, err:
        usage()
        print str(err), "\n"
        exit(2)

    print " "

    opts = {}
    for o, a in options: 
        opts[o] = a
    if "-v" in opts:
        DEBUG = True

    if '-R' in opts:
        recursive = True
        if DEBUG: print "recursive : Yes"
    if '-m' in opts:
        try:
            mode = int(opts['-m'])
            if DEBUG: print "%d pose(s) will be selected" % mode
        except:
            print "ERROR in specifying the poses number:", opts['-m']
            usage()
            exit(1)
    if '-p' in opts:
        pattern = opts['-p']
        #print "PATTERN", pattern
        if DEBUG: print "File pattern:", pattern
    if '-f' in opts:
        input_files = opts['-f'].split(",")
        if DEBUG: print "[ file mode:", input_files,"]"
    elif '-d' in opts:
        dir_root = opts['-d']
        if DEBUG: print "[ directory mode: %s ]" % dir_root
        if recursive:
            print "- scanning sub-directories of [ %s ]..." % dir_root
        else:
            print "- scanning directory [ %s ]" % dir_root
        #print recursive, pattern
        input_files = pathToList(dir_root, recursive = recursive, pattern = pattern)
        print "- %d results found" % (len(input_files))
        #print input

    else:
        print "ERROR: Missing input source: specify either file or directories"
        usage()
        exit(2)
    if not input_files:
        print "No input files found."
        exit(0)
    if '-r' in opts: # TODO caching could occur here?
        receptor = opts['-r']
        if not os.path.isfile(opts['-r']):
            print "ERROR! the specified filename is not accessible:", receptor 
            exit(1)
    else:
        print "No receptor specified (interactions will not be calculated)"
        doInteractions = False
    if '-i' in opts:
        doInteractions = False
    if '-o' in opts:
        suffix = opts['-o']
        if DEBUG: print "changed suffix:", suffix        
    if '-l' in opts:
        logfile = opts['-l']
        try:
            if os.path.exists(logfile):
                log = open(logfile, 'a')
            else:
                log = open(logfile, 'w')
                log.write(header)
        except:
            print "ERROR opening the log file:", logfile
            exit(1)

    if '-x' in opts:
        try:
            hbtol = float(opts['-x'])
        except:
            print "ERROR setting hydrogen bond extra tolerance."
            usage()
            exit(1)
            
    print "---------------------------------------------------------"


    ######################### start here
    if DEBUG: 
        print "[ receptor : %s ]" % receptor
        print "[ %d input files to process ]" % len(input_files)
        if doInteractions: print "[ interactions are calculated ]"
        else: print "[ interactions are NOT calculated ]"
    c = 0
    for l in input_files:
        try:
            c+=1
            if DEBUG: print "processing", l
            lig =  AutoDockVinaVsResult(input_files = l, 
                    mode = mode,
                    receptor = receptor, 
                    recname = None, 
                    auto = True, 
                    doInteractions = doInteractions,
                    hbtol = hbtol)
            pdbqt = lig.generatePDBQTplus()
            path_root = os.path.dirname(l) # the filename will be used to extract the root dir for the output file!
            if not path_root:
                path_root = os.getcwd()
            output_filename = path_root+os.sep+lig.ligName+suffix+'.pdbqt' 
            fp = open(output_filename, 'w')
            fp.write(pdbqt)
            fp.close()

            out_line = "%s\t%2.3f\t%2.3f\t%d\t%s\n" % ( lig.ligName, lig.results[0]['energy'],
                                               lig.results[0]['leff'], len(lig.poses), l)

            if '-l' in opts:
                print "\r- processing ... [ %d | %d  ]                     " % (c, len(input_files)),
                log.write(out_line)           
            else:
                print out_line[:-1]
            #if DEBUG and doInteractions:
            #    doDebug(l)
        except:
            print "ERROR: problems processing input :"
            print "file :", l
            print "error:", exc_info()[1]
    if '-l' in opts:
        log.close()
    print "\rDone                                            "

    exit(0)

