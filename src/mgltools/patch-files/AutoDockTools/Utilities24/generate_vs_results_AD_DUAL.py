#!/usr/bin/python

#############################################################################
#
# Author: Stefano FORLI
#
# Copyright: Stefano Forli, TSRI 2011
#
# v.0.2
#############################################################################


# specific helper functions
def findLigandsInDlg(dlg_list, checkSuccess = True):
    """

    """
    ligands = {}
    problematic = []
    for f in dlg_list:
        try:
            lines = getLines(f)
            for l in lines:
                if l.startswith("DPF> move"):
                    l = l.split("DPF> move", 1)[1].split(".pdbqt")[0]
                    if checkSuccess and not "Successful" in lines[-5]:
                        raise
                    if not l in ligands:
                        ligands[l] = [f]
                    else:
                        ligands[l].append(f)
                    break
        except:
            if DEBUG: print "[debug] problem reading file:",f
            error = sys.exc_info()[1]
            problematic.append([f, error])
    return ligands, problematic

def doDebug(lig):
    print lig
    if doInteractions:
        buffer = []
        try:
        #if 1:
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
        # else: 
            print "ERROR! writing the debug files"
            exit(2)


def usage():
    print """

INFO     Process AutoDock virtual screening results generating the following output:
             - PDBQT+ files  coordinate files containing interactions and properties
             - log file      summary of the docking results sorted by energy


USAGE    %s  [-f filename|-d directory] -r receptor [ OPTIONS ]

            -f filename     DLG file(s) to process. Multiple files can be separated by commas:
                            "ligand1_protein.dlg,ligand2_protein.dlg,..."
            -d directory    directory containing DLG files to process
            -r receptor     specify receptor filename (PDBQT)

OPTIONS
            -i              disable ligand-target interactions calculation 
            -o suffix       override default suffix "_VS" after ligand name in the output file
            -l filename     write ligand results summary to comma separated values (CSV) filename
                               - if the file does not exist, it will be created
                               - if the file exists, results ligand results will be append
            -m [all,le,lc]  specify which pose (LE, LC) properties are reported the log  [ default: all ]
            -t X.XX         set the RMSD tolerance for clustering [ default 2.00 Angstroms ]
            -R              scan all subdirectories in the specified directory
            -A              disable check for DLG completion before importing
            -h X.XX         set the extra tolerance for hb distance

            -x filename     calculate and print RMSD by using this reference PDBQT structure
                            (hydrogens ignored)

            -v              enable verbose mode (DEBUG)


EXAMPLE
            %s -d . -r protein.pdbqt -t 1.0 

            Recluster all DLGs in the current directory with 1.0A RMSD tolerance and use "protein.pdbqt"
            coordinates to calculate ligand-target interactions.
""" % (argv[0], argv[0])

#            -R              scan specified directory and all sub-directories for DLG(s)
#                            ( ignored with -f)

#===============================================
if __name__ == '__main__':
    #from glob import glob  
    import os
    from sys import argv, exc_info
    import getopt
    from AutoDockTools.VsResultsGenerator import *
    from AutoDockTools.HelperFunctionsN3P import *
    #from AutoDockTools.WaterProcessing import *


    # Defaults
    DEBUG = False
    doInteractions = True
    recursive = False
    rmsTol = 2.0
    #default_mode = 1 # 0:LE, 1:LC
    requested_mode = 'all' # by default, both LE and LC are printed in the log file
    pattern = "*.dlg"

    # clean up and prepare problematic log
    probl_log = "problematic_files.log"
    try: os.remove(probl_log)
    except: pass
    FOUND_PROBLEMS=0

    #header = "#name,energy,l_efficiency,clust_size,clust_size_percent,total_poses,filename,is_hydrated\n"
    header = "#name\tpose\tenergy\tl_efficiency\tclust_size\tclust_size_percent\ttotal_poses\tis_hydrated\tfilename"
    suffix = ".VS"
    receptor = None
    water_map = None
    checkSuccess = True
    fullentropy=False
    hbtol = 0.0
    try:
        options, var = getopt.getopt(argv[1:], 'f:d:r:o:t:ivm:Al:RAh:x:W:E')
        # -f  str       dlg single file
        # -d  str       directory
        # -r  str       receptor file
        # -R  (none)    recursive (find all dlgs in a given path)
        # -o  str       output filename (override defult "ligname_VS.pdbqt")
        # -O  dir       save the output PDBQT+ file in this path? TODO TODO TODO
        # -m  str       set the mode for which energy data must be printed ("le", "lc")
        # -A            disable DLG success check
        # -l  logfile   add the logfile info
        # -h  float     add the extra tolerance to the hb distance
        # -x  filename  calculate RMSD with xray structure
        # 
        # -t  int       RMS tolerance
        # -i            ignore interactions
        # -v  (none)    verbose (DEBUG MODE)
        # -W  mapfile   water map file
        # -E            full entropy treatment (penalty for conserved waters) XXX REMOVE, PROBABLY USELESS

        # TODO very important!!!
        # TODO support the caching of the receptor? Why not?....

    except getopt.GetoptError, err:
        usage()
        print str(err), "\n"
        exit(2)
    opts = {}

    print " " 
    # parse options
    for o, a in options: 
        opts[o] = a
    if "-v" in opts:
        DEBUG = True
    if '-r' in opts:
        receptor = opts['-r']
        if not os.path.isfile(opts['-r']):
            print "ERROR! the specified filename is not accessible:", receptor 
            exit(1)
        print "- caching receptor structure '%s'..." % receptor,
        try:
            receptor=getCoords(getLines(receptor))
            print "[ Done ]"
        except:
            print "ERROR loading receptor structure [%s]" % receptor
            exit(2)
    else:
        print "No receptor specified (interactions will not be calculated)"
        doInteractions = False
    if '-i' in opts:
        doInteractions = False
        print "Interactions will be *not* calculated"


    if '-A' in opts:
        checkSuccess = False
        if DEBUG: print "Disable DLG completion check (\"Successful\")"
    if '-R' in opts:
        recursive = True
    if '-m' in opts:
        if (opts['-m'] == 'le'):
            requested_mode = [0]
        elif (opts['-m'] == 'lc'):
            requested_mode = [1]
        else:
            print "ERROR in specifying the pose mode (le/lc)", opts['-m']
            usage()
            exit(1)
    if '-f' in opts:
        input_file = opts['-f'].split(",")
        if DEBUG:
            print "[ file mode:", input_file,"]"
        input_file, problematic = findLigandsInDlg(input_file, checkSuccess = checkSuccess)
        print " [ %d ligands, %d skipped/problematic ]" % (len(input_file), len(problematic))
    elif '-d' in opts:
        dir_root = os.path.abspath(opts['-d'] )
        #opts['-d']
        if DEBUG: print "[ directory mode: %s ]" % dir_root
        if recursive:
            print "- scanning sub-directories in '%s'..." % dir_root,
        else:
            print "- scanning directory '%s'" % dir_root,
        if not os.path.exists(dir_root):
            print "ERROR: directory doesn't exist"
            usage()
            exit(2)
        print "[ Done ]"
        dlg_list = pathToList(dir_root, recursive = recursive, pattern = pattern)
        print "- %d DLG found" % (len(dlg_list)),
        input_file, problematic = findLigandsInDlg(dlg_list, checkSuccess = checkSuccess)
        print " [ %d ligands, %d skipped/problematic ]" % (len(input_file), len(problematic))
    else:
        print "ERROR: Missing input source: specify either file or directories"
        usage()
        exit(2)
    if not input_file:
        print "No input files found." 
        exit(0)
    if '-t' in opts:
        try:
            rmsTol = float(opts['-t'])
        except:
            print "ERROR setting RMSD tolerance to '%s'" % opts['-t']
            usage()
            exit(2)
    if '-o' in opts:
        suffix = opts['-o']
        if DEBUG: print "changed suffix:", suffix
    if '-l' in opts:
        logfile = opts['-l']
        print "- log file : %s " % logfile,
        try:
            if os.path.exists(logfile):
                log = open(logfile, 'a')
                print "[ appending ]"
            else:
                log = open(logfile, 'w')
                log.write(header+'\n')
                print "[ creating new ]"
        except:
            error = sys.exc_info()[1]
            print "ERROR opening the log file:", logfile, 
            print error
            exit(1)

    if '-h' in opts:
        try:
            hbtol = float(opts['-h'])
            if DEBUG: print "[ Extra hb tolerance : %1.3f ]" % (hbtol)
        except:
            print "ERROR setting extra hb tolerance."
            usage()
            exit(1)

    if '-W' in opts:
        try:
        #if True:
            data = getLines(opts['-W'])
            water_map = map2array(data)

        except:
            print "ERROR when loading the water map [ %s ]" % opts['-W']
            exit(1)
    
    if '-E' in opts:
        fullentropy=1
        print '- *FULL* entropy treatment' #XXX REMOVE probably useless

    # end parsing options


    # Problems pre-scanning report:
    if problematic:
        FOUND_PROBLEMS=1
        print "[ => problematic scanned files are saved in log file : %s ]" % probl_log
        fp = open(probl_log, 'a')
        #fp.write('== Problems scanning input files ==\n')
        for p in problematic:
            line = "SCAN_ERROR>%s : %s\n" % (p[0],p[1])
            fp.write(line)
        fp.write('===================================\n')
        fp.close()

    ######################### start here
    if DEBUG: 
        print "[ receptor : %s ]" % receptor
        print "[ %d input files to process ]" % len(input_file)
        if doInteractions:
            print "[ interactions are calculated ]"
        else:
            print "[ interactions are NOT calculated ]"
        print "[ RMSD tolerance : %2.1f ]" % rmsTol
    c = 0

    # TODO show the header of the data if '-l' not specified?
    print "---------------------------------------------------------"

    result_sort = []
    for l in input_file.keys():
        #mode = requested_mode
        #try: XXX CHANGE THIS IN THE FINAL VERSION!!
        if 1:
            c += 1
            if DEBUG: print "processing", input_file[l]
            lig =  AutoDockVsResult(input_files = input_file[l], 
                    rmsTol = rmsTol,
                    receptor = receptor, 
                    recname = None, 
                    auto = True, 
                    doInteractions = doInteractions,
                    water_map = water_map,
                    hbtol=hbtol,
                    fullentropy=fullentropy)

            # save PDBQT+
            pdbqt = lig.generatePDBQTplus()
            # output file naming 
            path_root = os.path.dirname(input_file[l][0]) # the first DLG of the ligand is used to get the path name
            if not path_root:
                path_root = os.getcwd()
            output_filename = path_root+os.sep+lig.ligName+suffix+'.pdbqt'
            fp = open(output_filename, 'w')
            fp.write(pdbqt)
            fp.close()


            # TODO think about: what if user wants to print *both* ?
            hydrated = "\t"
            if lig.isHydrated:
                hydrated = "\tHYDRATED"

            # XXX req mode here
            if requested_mode == 'all':
                pool = range(len(lig.results))
            else:
                if requested_mode+1 > len(lig.results):
                    pool = [0]
                else:
                    pool = requested_mode


            for p in pool:
                if len(lig.results) == 1 and mode == 1: # LC requested but only one result (LE == LC)
                #    mode = 0
                #    if DEBUG: print "[ only one result found, switching to LE mode ]"
                #if len(lig.results) > 0:
                pose = lig.results[p]
            #else:
            #    print "No poses found!"
            #    raise
            

            #out_line = "%s,%2.3f,%2.3f,%d,%2.1f,%d,%s%s\n" % ( lig.ligName, pose['energy'],
                out_line = "%s\t%d\t%2.3f\t%2.3f\t%d\t%2.1f\t%d\t%s%s\n" % ( lig.ligName, p+1, pose['energy'],
                                                           pose['leff'], pose['csize'],
                                                           pose['cpercent'], len(lig.poses),
                                                           hydrated,output_filename)
                result_sort.append( (lig.ligName, p+1, pose['energy'],
                                                           pose['leff'], pose['csize'],
                                                           pose['cpercent'], len(lig.poses),
                                                           hydrated,output_filename) )
                if '-l' in opts:
                    log.write(out_line)
                    if not c % 10:
                        print "\r- processing ... [ %d | %d  ]                     " % (c, len(input_file)),
                else:
                    print out_line[:-1]
                if DEBUG:
                    doDebug(lig)

                if '-x' in opts:
                    print "RMSD_reference [%s]" % opts['-x']
                    file_list = []
                    rms,dist_pool,poses = lig.calcReferenceRms(opts['-x'], debug=True)
                    print len(dist_pool)
                    for p in range(len(dist_pool)):
                        xxx = []
                        tot = 0
                        for c in range(len(dist_pool[p])):
                            t = "%d %2.3f" % (c, dist_pool[p][c])
                            xxx.append(t)
                            tot += dist_pool[p][c]
                        name = "pose_rmsd_pool_%d" % (p) # , tot)
                        log = "%2.3f\t %s" % (tot, name)
                        file_list.append(log)
                        writeList(name+".dat", xxx, addNewLine=True)
                        writeList(name+".pdb", poses[p]['text'])
                        
                    writeList('log_list.log', file_list, addNewLine=True)
                    print "RMSD_LE> %2.3f A" % rms[0]
                    if len(rms) >1:
                        print "RMSD_LC> %2.3f A" % rms[1]


        #except:
        else:
            # debugging printing this?
            if DEBUG:
                print "\nWARNING: problems processing one or more of the input files"
                print "    File(s)    : %s "% (input_file[l])
                print "    Error [%s] : %s "% (exc_info()[0], exc_info()[1])
            FOUND_PROBLEMS=1
            fp = open(probl_log,'a')
            txt = "PROCESS_ERROR>%s : [%s] %s" % (input_file[l], exc_info()[0], exc_info()[1])
            fp.write(txt)
            fp.close()
            # TODO add the problematic file here...

    if '-l' in opts:
        log.close()

        # XXX overwrite the temp file and write the sorted one
        result_sort = sorted(result_sort, key=lambda x: x[2])
        fp = open(logfile, 'w')
        fp.write(header+'\n')
        for r in result_sort:
            out_line = "%s\t%d\t%2.3f\t%2.3f\t%d\t%2.1f\t%d\t%s%s\n" % ( r )
            fp.write(out_line)
        fp.close()

    if FOUND_PROBLEMS:
        print "\nWARNING: Some files raised errors when processing."
        print "         and have been saved in:\n"
        print "         %s" % probl_log
    print "\rDone                                            "
    exit(0)
