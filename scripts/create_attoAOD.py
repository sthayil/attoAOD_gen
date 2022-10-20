#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import glob, argparse, socket
from datetime import datetime

#===================== define input/output here ===========================
OutputFile = ""      #if set to "", then OutputFile = <FileListName> + _plots.root
MaxEvents = 0        #set to 0 to run all events in each file
#MaxFiles = -1        #set to -1 to run all files in file list
#MaxFiles = 500

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] outputDir inputFiles")

    parser.add_option("-m", "--maxFiles", dest="MaxFiles", type=int, default=500,
                      help="Max #files to run over from filelist")
    parser.add_option("--batch", dest="Batch", type=int, default=None,
                      help="Batch number of this run")
    parser.add_option("-l", "--lepton", dest="lepton", type="string", default='mu',
                      help="Lepton to run over (el/mu)")
    parser.add_option("-f", "--filelist", dest="FileList", type="string", default=None,
                      help="Filelist to run over")

    parser.add_option("-s", "--postfix", dest="postfix", type="string", default=None,
                      help="Postfix which will be appended to the file name (default: _Friend for friends, _Skim for skims)")
    parser.add_option("-J", "--json", dest="json", type="string",
                      default=None, help="Select events using this JSON file")
    parser.add_option("-c", "--cut", dest="cut", type="string",
                      default=None, help="Cut string")
    parser.add_option("-b", "--branch-selection", dest="branchsel",
                      type="string", default=None, help="Branch selection")
    parser.add_option("--bi", "--branch-selection-input", dest="branchsel_in",
                      type="string", default=None, help="Branch selection input")
    parser.add_option("--bo", "--branch-selection-output", dest="branchsel_out",
                      type="string", default=None, help="Branch selection output")
    parser.add_option("--friend", dest="friend", action="store_true", default=False,
                      help="Produce friend trees in output (current default is to produce full trees)")
    parser.add_option("--full", dest="friend", action="store_false", default=False,
                      help="Produce full trees in output (this is the current default)")
    parser.add_option("--noout", dest="noOut", action="store_true",
                      default=False, help="Do not produce output, just run modules")
    parser.add_option("-P", "--prefetch", dest="prefetch", action="store_true", default=False,
                      help="Prefetch input files locally instead of accessing them via xrootd")
    parser.add_option("--long-term-cache", dest="longTermCache", action="store_true", default=False,
                      help="Keep prefetched files across runs instead of deleting them at the end")
    parser.add_option("-N", "--max-entries", dest="maxEntries", type="long", default=None,
                      help="Maximum number of entries to process from any single given input tree")
    parser.add_option("--first-entry", dest="firstEntry", type="long", default=0,
                      help="First entry to process in the three (to be used together with --max-entries)")
    parser.add_option("--justcount", dest="justcount", default=False,
                      action="store_true", help="Just report the number of selected events")
    parser.add_option("-I", "--import", dest="imports", type="string", default=[], action="append",
                      nargs=2, help="Import modules (python package, comma-separated list of ")
    parser.add_option("-z", "--compression", dest="compression", type="string",
                      default=("LZMA:9"), help="Compression: none, or (algo):(level) ")

    (options, args) = parser.parse_args()

    if options.friend:
        if options.cut or options.json:
            raise RuntimeError(
                "Can't apply JSON or cut selection when producing friends")

    modules = []
    for mod, names in options.imports:
        import_module(mod)
        obj = sys.modules[mod]
        selnames = names.split(",")
        mods = dir(obj)
        for name in selnames:
            if name in mods:
                print("Loading %s from %s " % (name, mod))
                if type(getattr(obj, name)) == list:
                    for mod in getattr(obj, name):
                        modules.append(mod())
                else:
                    modules.append(getattr(obj, name)())
    if options.noOut:
        if len(modules) == 0:
            raise RuntimeError(
                "Running with --noout and no modules does nothing!")
    if options.branchsel != None:
        options.branchsel_in = options.branchsel
        options.branchsel_out = options.branchsel

    def read_file_list(FileList, MaxFiles, Batch):
        f=open(FileList, "r")
        InputFiles = f.readlines()
        f.close()

        nInputFiles = len(InputFiles)

        startfile = int(Batch)*int(MaxFiles)
        endfile = startfile + int(MaxFiles)
        if nInputFiles<startfile:
            raise RuntimeError(
                "Batch number is too high!")
        if nInputFiles<endfile: endfile=nInputFiles

        # if MaxFiles > nInputFiles:
        #     print "MaxFiles", MaxFiles, "> nInputFiles", nInputFiles
        #     quit()

        # nOutputFiles = MaxFiles
        # if MaxFiles == -1: nOutputFiles = nInputFiles

        hostname = socket.gethostname()
        if ".fnal.gov" in hostname:
            for i in range(startfile,endfile): InputFiles[i] = 'root://cmseos.fnal.gov/'+InputFiles[i].strip()
        elif "hexcms" in hostname:
            for i in range(startfile,endfile): InputFiles[i] = InputFiles[i].strip(' \n')
        return InputFiles[startfile:endfile]

    FileList = options.FileList
    Batch = options.Batch
    if OutputFile == "":
        OutputFile = FileList.split("/")[-1]
        OutputFile = OutputFile.split(".")[0]
        OutputFile = str(Batch)+"_"+options.lepton+"_"+ OutputFile + ".root"

    print read_file_list(FileList, options.MaxFiles, Batch)

    p = PostProcessor(".", read_file_list(FileList, options.MaxFiles, Batch),
                      cut=options.cut,
                      branchsel=options.branchsel_in,
                      modules=modules,
                      compression=options.compression,
                      friend=options.friend,
                      haddFileName=OutputFile,
                      #postfix=options.postfix,
                      jsonInput=options.json,
                      noOut=options.noOut,
                      justcount=options.justcount,
                      prefetch=options.prefetch,
                      longTermCache=options.longTermCache,
                      maxEntries=options.maxEntries,
                      firstEntry=options.firstEntry,
                      outputbranchsel=options.branchsel_out)
    p.run()
