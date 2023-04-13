# use condor to make attoAODs from nanoAODs
import os, glob, sys, argparse, socket
from datetime import datetime
options = argparse.ArgumentParser(description="Sets up a run to generate attoAODs from nanoAODs")
options.add_argument("-l", "--lepton",           required=True, help="lepton to be used for analysis", choices=['el','mu'])
options.add_argument("-y", "--year",             required=True, help="year from run2", choices=['2016','2017','2018'])
options.add_argument("-d", "--dataset",          required=True, help="dataset", choices=['dyjetstoll','wjetstolnu','ttjets','egammaA','egammaB','egammaC','egammaD','singlemuonA','singlemuonB','singlemuonC','singlemuonD','sig_M1000','sig_M500','sig_M250'])
options.add_argument("-n", "--numFiles",         nargs='?',     help="#files from filelist in each job (defaults to 500)", const=500, type=int, default=500)
options.add_argument("-o", "--outputDirectory",  nargs='?',     help="Output base directory filepath for jobs. Should be an EOS area. (Defaults to /store/user/lpcrutgers/sthayil/pseudoaxions/atto_passTrigger)", const='/store/user/lpcrutgers/sthayil/pseudoaxions/atto_passTrigger', type=str, default='/store/user/lpcrutgers/sthayil/pseudoaxions/atto_passTrigger')
options.add_argument("-m", "--mode",             nargs='?',     help="Run mode", const='normal', type=str, default='normal', choices=['normal','hadd'])
options.add_argument("-hd", "--haddDir",         nargs='?',     help="Directory to put hadded files in", const='hadded', type=str, default='hadded')
ops = options.parse_args()

#hadd mode--------------------------------------------------------------------------------
if ops.mode=='hadd':
    if not os.path.isdir(ops.haddDir): 
        print "Making job directory "+ops.haddDir
        os.system('mkdir '+ops.haddDir)
    hostname = socket.gethostname()
    jobname = ops.lepton+'_'+ops.dataset+'_'+ops.year
    if ".fnal.gov" in hostname:
        if (ops.outputDirectory).startswith("/store/user/"):
            os.system('eos root://cmseos.fnal.gov ls '+ops.outputDirectory+'/'+jobname+' > haddfilelist')
            with open('haddfilelist', 'r+') as in_file:
                rootFiles=[]
                totEntrFiles=[]
                for count, line in enumerate(in_file):
                    if ".root" in line: rootFiles.append(line.rstrip('\n'))
                    elif ".txt" in line: totEntrFiles.append(line.rstrip('\n'))
                if len(rootFiles)!=len(totEntrFiles): 
                    print "#root files != #totentries files; check what's going on"
                    exit()

                os.system('rm haddfilelist')
                os.chdir(jobname)

                command='python ../scripts/haddnano.py '+jobname+'.root '
                for rootFile in rootFiles: 
                    os.system('xrdcp -s root://cmsxrootd.fnal.gov/'+ops.outputDirectory+'/'+jobname+'/'+rootFile+' .')
                    command+=rootFile+' '
                os.system(command)
                os.system('mv '+jobname+'.root ../'+ops.haddDir)
                for rootFile in rootFiles: os.system('rm '+rootFile)
                tottotEntries=0
                for totEntrFile in totEntrFiles: 
                    os.system('xrdcp -s root://cmsxrootd.fnal.gov/'+ops.outputDirectory+'/'+jobname+'/'+totEntrFile+' .')
                    with open(totEntrFile, 'r') as in_file: 
                        dataline=in_file.readline()
                        tottotEntries+=int(dataline)
                    os.system('rm '+totEntrFile)
                print "\n\n"+str(tottotEntries)
                with open('../'+ops.haddDir+'/'+jobname+'_totentries.txt', 'w') as entriesfile:
                    entriesfile.write(str(tottotEntries))
        else:
            print "ERROR: Output directory is not an EOS path starting in /store/user/"
            exit()
    elif "hexcms" in hostname: 
        if os.path.isdir(ops.outputDirectory+'/'+jobname):
            os.system('ls '+ops.outputDirectory+'/'+jobname+' > haddfilelist')
            with open('haddfilelist', 'r+') as in_file:
                rootFiles=[]
                totEntrFiles=[]
                for count, line in enumerate(in_file):
                    if ".root" in line: rootFiles.append(line.rstrip('\n'))
                    elif ".txt" in line: totEntrFiles.append(line.rstrip('\n'))
                if len(rootFiles)!=len(totEntrFiles): 
                    print "#root files != #totentries files; check what's going on"
                    exit()

                os.system('rm haddfilelist')
                os.chdir(jobname)

                command='python ../scripts/haddnano.py '+jobname+'.root '
                for rootFile in rootFiles: 
                    #os.system('xrdcp -s root://cmsxrootd.fnal.gov/'+ops.outputDirectory+'/'+jobname+'/'+rootFile+' .')
                    command+=ops.outputDirectory+'/'+jobname+'/'+rootFile+' '
                os.system(command)
                os.system('mv '+jobname+'.root ../'+ops.haddDir)
                #for rootFile in rootFiles: os.system('rm '+rootFile)
                tottotEntries=0
                for totEntrFile in totEntrFiles: 
                    #os.system('xrdcp -s root://cmsxrootd.fnal.gov/'+ops.outputDirectory+'/'+jobname+'/'+totEntrFile+' .')
                    with open(ops.outputDirectory+'/'+jobname+'/'+totEntrFile, 'r') as in_file: 
                        dataline=in_file.readline()
                        tottotEntries+=int(dataline)
                    #os.system('rm '+totEntrFile)
                print "\n\n"+str(tottotEntries)
                with open('../'+ops.haddDir+'/'+jobname+'_totentries.txt', 'w') as entriesfile:
                    entriesfile.write(str(tottotEntries))
        else:
            print "ERROR: Output directory: "+ops.outputDirectory+'/'+jobname+" is not a valid path"
    exit()

#directory with log files--------------------------------------------------------------------
jobname = ops.lepton+'_'+ops.dataset+'_'+ops.year
if not os.path.isdir(jobname): 
    print "Making job directory "+jobname
    os.system('mkdir '+jobname)

#directory (eos area) with attoAODs----------------------------------------------------------
hostname = socket.gethostname()
prefix=""
if ".fnal.gov" in hostname: 
    if (ops.outputDirectory).startswith("/store/user/"): 
        os.system('eos root://cmseos.fnal.gov mkdir -p '+ops.outputDirectory+'/'+jobname)
    # elif (ops.outputDirectory).startswith("/eos/uscms/store/user/") : #this wont work with the prefix bit; fix---------------
    #     os.system('eos root://cmseos.fnal.gov mkdir -p '+ops.outputDirectory+'/'+jobname)
    else:
        #print "ERROR: for output directory, specify an EOS path starting in /store/user/ or /eos/uscms/store/user/"
        print "ERROR: for output directory, specify an EOS path starting in /store/user/"
        exit()
    prefix="root://cmseos.fnal.gov/"
elif "hexcms" in hostname: 
    if os.path.isdir(ops.outputDirectory):
        if not ( (ops.outputDirectory).startswith("/") ): prefix=os.getcwd()+"/" #for hex, don't care about eos area, can save ouput locally
    else:
        print "provided output directory does not exist, trying to create..."
        if not ( (ops.outputDirectory).startswith("/") ): prefix=os.getcwd()+"/"
        os.system('mkdir -p '+ops.outputDirectory)
        print "created directory: "+prefix+ops.outputDirectory
    if os.path.isdir(ops.outputDirectory+'/'+jobname): 
        print "WARNING: directory named "+ops.outputDirectory+'/'+jobname+" already exists, files will be overwritten"
    else: os.system('mkdir -p '+ops.outputDirectory+'/'+jobname)
outputDir=prefix+ops.outputDirectory+'/'+jobname

#check accessibility of filelists-----------------------------------------------------------
filelist = 'filelists/'+ops.dataset+'_'+ops.year+'.txt' 
if not os.path.exists(filelist): 
    print "Check that your filelist exists (format: filelists/dataset_year.txt)"
    exit
f = open(filelist, "r")
flines = f.readlines()
f.close()

if ".fnal.gov" in hostname and not (flines[0].startswith('/store/user/')): print "Filelist items don't start in /store/user/; first item is: "+flines[0]
elif "hexcms" in hostname and not (flines[0].startswith('/cms/')): print "Filelist items don't start in /cms/; first item is: "+flines[0]

#print #files/job, cd to jobdir, calculate num jobs------------------------------------------
print "Each job will run over up to "+str(ops.numFiles)+" NanoAODs"

os.chdir(jobname)
if not os.path.isdir('logs'): os.system('mkdir logs')
if not os.path.isdir('jdl_files'): os.system('mkdir jdl_files')

numjobs = len(flines)//ops.numFiles + (len(flines)%ops.numFiles >0)

#modify jdl---------------------------------------------------------------------------------
mytime=datetime.now()
mytimestr=mytime.strftime("%y%m%d%H")
with open('../condorsubmit_nanotoatto.jdl', 'r') as file :
    filedata = file.read()
filedata = filedata.replace('00000000', ops.lepton) #ops.lepton ops.dataset ops.year outputDir ops.numFiles numjobs
filedata = filedata.replace('11111111', ops.dataset)
filedata = filedata.replace('22222222', ops.year)
filedata = filedata.replace('33333333', outputDir)
filedata = filedata.replace('44444444', str(ops.numFiles))
filedata = filedata.replace('55555555', str(numjobs))
with open('condorsubmit_lhetominiaod_'+mytimestr+'.jdl', 'w') as file:
    file.write(filedata)
    
#submit jobs-------------------------------------------------------------------------------
os.system('cp ../condorsubmit_nanotoatto.sh .')
print "Submitting jobs with condor_submit condorsubmit_nanotoatto_"+mytimestr+".jdl ..."
os.system('condor_submit condorsubmit_lhetominiaod_'+mytimestr+'.jdl')
os.system('mv condorsubmit_lhetominiaod_'+mytimestr+'.jdl jdl_files/')

#totalEvents mode-------------------------------------------------------------------------
