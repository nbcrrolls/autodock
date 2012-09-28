#!/usr/bin/python

#############################################################################
#
# Author: Stefano FORLI
#
# Copyright: Stefano Forli, TSRI 2011
#
# v.0.4
#############################################################################


# specific helper functions
def findLigandsInDlg(dlg_list, checkSuccess = True):
    """

    """
    from AutoDockTools.HelperFunctionsN3P import QuickStop
    ligands = {}
    problematic = []
    error = None
    for f in dlg_list:
        try:
            lines = getLines(f)
            for l in lines:
                if l.startswith("DPF> move"):
                    l = l.split("DPF> move", 1)[1].split(".pdbqt")[0]
                    if checkSuccess and not "Successful" in lines[-5]:
                        #raise "DOCKING" # Exception 'Docking not successful' # , f) # QuickStop
                        error = "Docking not successful"
                        break
                    if not l in ligands:
                        ligands[l] = [f]
                    else:
                        ligands[l].append(f)
                    break
        except:
            if DEBUG: print "[debug] problem reading file:",f
            error = sys.exc_info()[1]
        if not error == None:
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

            -q              suppress printing progress status
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
    import os, getopt
    import  time # DEBUG ONLY?
    from sys import argv, exc_info

    from AutoDockTools.VsResultsGenerator import *
    from AutoDockTools.HelperFunctionsN3P import *
    from AutoDockTools.WaterProcessing import *

    version = "0.4"

    # Defaults
    DEBUG = False
    QUIET = False
    doInteractions = True
    recursive = False
    rmsTol = 2.0
    # default_mode = 1 # 0:LE, 1:LC
    requested_mode = [0,1]
    pattern = "*.dlg"
    #header = "#name,energy,l_efficiency,clust_size,clust_size_percent,total_poses,filename,is_hydrated\n"
    header = "#name\tpose\tenergy\tl_efficiency\tclust_size\tclust_size_percent\ttotal_poses\tfilename\tis_hydrated"
    suffix = ".VS"
    receptor = None
    water_map = None
    checkSuccess = True
    fullentropy=False
    hbtol = 0.0

    print "\n\t * VS_GENERATOR V.%s" % version

    try:
        options, var = getopt.getopt(argv[1:], 'f:d:r:o:t:ivm:Al:RAh:x:W:Eq')
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

    if '-q' in opts:
        QUIET = True

    if '-r' in opts: # TODO caching could occur here?
        receptor = opts['-r']
        if not os.path.isfile(opts['-r']):
            print "ERROR! the specified filename is not accessible:", receptor 
            exit(1)
        print "- receptor file : %s " % receptor
    else:
        print "No receptor specified (interactions will not be calculated)"
        doInteractions = False
    if '-i' in opts:
        doInteractions = False
        print "Interactions will be *not* calculated"


    if '-A' in opts:
        checkSuccess = False
        if DEBUG: print "Disable DLG check (\"Successful\")"

    if '-R' in opts:
        recursive = True

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

    if '-m' in opts:
        if (opts['-m'] == 'le'):
            requested_mode = [0]
        elif (opts['-m'] == 'lc'):
            requested_mode = [1]
        elif opts['-m'] == 'all':
            requested_mode = [0,1]
        else:
            print "ERROR in specifying the pose mode (le/lc)", opts['-m']
            usage()
            exit(1)
    else:
        if requested_mode == [0]:
            mode = ' 0\t[ lowest energy ]'
        elif requested_mode == [1]:
            mode = ' 1\t[ largest cluster ]'
        else:
            mode = ' 2\t[ lowest energy + largest cluster ]'
        print "- energy mode : %s " % mode

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


    # here are the time-consuming parts
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
            print "- scanning sub-directories of [ %s ]..." % dir_root
        else:
            print "- scanning directory [ %s ]" % dir_root
        dlg_list = pathToList(dir_root, recursive = recursive, pattern = pattern)
        print "- %d DLG found" % (len(dlg_list)),
        input_file, problematic = findLigandsInDlg(dlg_list, checkSuccess = checkSuccess)
        print " [ %d ligands, %d skipped/problematic ]" % (len(input_file), len(problematic))
    else:
        print "ERROR: Missing input source: specify either file or directories"
        usage()
        exit(2)
    # end parsing options


    # Problems pre-scanning report:
    if len(problematic):
        probl_log = "problematic_results.log"
        print "[ => problematic result files are saved in log file : %s ]" % probl_log
        fp = open(probl_log, 'w')
        for p in problematic:
            line = "%s => %s\n" % (p[0],p[1])
            fp.write(line)
        fp.close()

    if not input_file:
        print "No input files found." 
        exit(0)


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
    #print "---------------------------------------------------------"

    result_sort = []

    processor = AutoDockVsResult(rmsTol = rmsTol,
                    receptor = receptor, 
                    recname = None, 
                    auto = True, 
                    doInteractions = doInteractions,
                    water_map = water_map,
                    hbtol=hbtol,
                    fullentropy=fullentropy)

    if QUIET: print "Processing...",
    t0 = time.time() 
    for l in input_file.keys():
        #mode = requested_mode
        #print input_file[l]
        processor.setLigands(input_file[l])

        #try: XXX CHANGE THIS IN THE FINAL VERSION!!
        if 1:
            c += 1
            if DEBUG: print "processing", input_file[l]
            """
            lig = AutoDockVsResult(input_files = input_file[l], 
                    rmsTol = rmsTol,
                    receptor = receptor, 
                    recname = None, 
                    auto = True, 
                    doInteractions = doInteractions,
                    water_map = water_map,
                    hbtol=hbtol,
                    fullentropy=fullentropy)
            """
            processor.process()
            pdbqt = processor.generatePDBQTplus()
            # output file naming 
            # XXX make this an option? a dir where p+ will be stored (split subdirs?)
            path_root = os.path.dirname(input_file[l][0]) # the first DLG of the ligand is used to get the path name
            if not path_root:
                path_root = os.getcwd()
            output_filename = path_root+os.sep+processor.ligName+suffix+'.pdbqt'
            fp = open(output_filename, 'w')
            fp.write(pdbqt)
            fp.close()

            # XXX req mode here
            if len(processor.results) == 0:
                print "Warning: No poses found!"
                raise Exception

            if requested_mode == [0]:
                pool = [0]
            elif requested_mode == [0,1]:
                pool = range(len(processor.results))
            elif requested_mode == [1]:
                if len(processor.results) == 2:
                    pool = requested_mode
                elif len(processor.results) == 1:
                    pool = [0]

            hydrated = ""
            if processor.isHydrated:
                hydrated = "\tHYDRATED"

            for p in pool:
                pose = processor.results[p]
                out_line = "%s\t%d\t%2.3f\t%2.3f\t%d\t%2.1f\t%d\t%s%s\n" % ( processor.ligName, p+1, pose['energy'],
                           pose['leff'], pose['csize'], pose['cpercent'], len(processor.poses),
                           output_filename,hydrated)

                result_sort.append( (processor.ligName, p+1,  pose['energy'], pose['leff'], pose['csize'],
                            pose['cpercent'], len(processor.poses), output_filename,hydrated) )

                if '-l' in opts:
                    log.write(out_line)
                    if not QUIET: print "\r- processing ... [ %d | %d  ]                     " % (c, len(input_file)),
                else:
                    if not QUIET:
                        print out_line[:-1]
                if DEBUG:
                    doDebug(lig)

                if '-x' in opts:
                    print "RMSD_reference [%s]" % opts['-x']
                    file_list = []
                    rms,dist_pool,poses = processor.calcReferenceRms(opts['-x'], debug=True)
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
            print "\nWARNING: problems processing one or more of the input files"
            print "    File  : %s "% (input_file[l])
            print "    Error : %s "% (exc_info()[1])

            fp = open(probl_log,'a')
            txt = "PROCESS_ERROR>%s : [%s] %s" % (input_file[l], exc_info()[0], exc_info()[1])
            fp.write(txt)
            fp.close()
    print "ELAPSED TIME:", time.time() - t0

    if '-l' in opts:
        log.close()
        print "- sorting results ...",
        # overwrite log file with energy-sorted results
        written_data = hf.getLines(logfile, doStrip=1)
        log_data = [ x.split() for x in written_data if not x.startswith("#") ]
        result_sort = []
        for r in log_data:
            r[1] = int(r[1])
            r[2] = float(r[2])
            r[3] = float(r[3])
            r[4] = int(r[4])
            r[5] = float(r[5])
            r[6] = int(r[6])
            result_sort.append(r)
        result_sort = sorted(result_sort, key=lambda x: x[2])
        fp = open(logfile, 'w')
        fp.write(header+'\n')
        for r in result_sort:
            if len(r) < 9:
                r.append("")
            out_line = "%s\t%d\t%2.3f\t%2.3f\t%d\t%2.1f\t%d\t%s\t%s" % ( tuple(r) )
            fp.write(out_line+'\n')
        fp.close()

        print "[ DONE ]"

    print "\rDone                                            "
    exit(0)
