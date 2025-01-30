# use condor to make attoAODs from nanoAODs
import os, glob, sys, argparse, socket,re
from datetime import datetime
options = argparse.ArgumentParser(description="Sets up a run to generate attoAODs from nanoAODs")
options.add_argument("-l", "--lepton",                          help="lepton to be used for analysis", choices=['el','mu'],default='mu')
options.add_argument("-y", "--year",             required=True, help="year from run2", choices=['2016','2016APV','2017','2018'])
options.add_argument("-d", "--dataset",          required=True, help="dataset", choices=['dyjetstoll','wjetstolnu','ttjets','singlemuonA','singlemuonB','singlemuonB2','singlemuonC','singlemuonD','singlemuonE','singlemuonF','singlemuonF2','singlemuonG','singlemuonH','eta_M500','eta_M750','eta_M850','eta_M1000','eta_M1500','eta_M2000','eta_M2500','eta_M3000','eta_M4000','etaprime_M850','etaprime_M1000','etaprime_M1500','etaprime_M2000','etaprime_M2500','etaprime_M3000','etaprime_M4000'])
options.add_argument("-n", "--numFiles",         nargs='?',     help="#files from filelist in each job (defaults to 500)", const=500, type=int, default=500)
options.add_argument("-o", "--outputDirectory",  nargs='?',     help="Output base directory filepath for jobs. Should be an EOS area. (Defaults to /store/user/lpcrutgers/sthayil/pseudoaxions/atto)", const='/store/user/lpcrutgers/sthayil/pseudoaxions/atto', type=str, default='/store/user/lpcrutgers/sthayil/pseudoaxions/atto')
options.add_argument("-m", "--mode",             nargs='?',     help="Run mode", const='normal', type=str, default='normal', choices=['normal','hadd'])
options.add_argument("-hd", "--haddDir",         nargs='?',     help="Directory to put hadded files in", const='hadded', type=str, default='hadded')
#options.add_argument("-s",  "--schedd",          nargs='?',     help="Specific scheduler to submit jobs to", type=str, default=None)
ops = options.parse_args()

#hadd mode--------------------------------------------------------------------------------
if ops.mode=='hadd':
    hostname = socket.gethostname()
    jobname = ops.lepton+'_'+ops.dataset+'_'+ops.year
    if ".fnal.gov" in hostname:
        # if not os.path.isdir(ops.haddDir): 
        #     print ("Making job directory "+ops.haddDir)
        #     os.system('mkdir '+ops.haddDir)

        if (ops.outputDirectory).startswith("/store/user/"):
            os.system('eos root://cmseos.fnal.gov ls '+ops.outputDirectory+'/'+jobname+' > haddfilelist')
            with open('haddfilelist', 'r+') as in_file:
                rootFiles=[]
                totEntrFiles=[]
                for count, line in enumerate(in_file):
                    if ".root" in line: rootFiles.append(line.rstrip('\n'))
                    elif ".txt" in line: totEntrFiles.append(line.rstrip('\n'))
                if len(rootFiles)!=len(totEntrFiles): 
                    print ("#root files != #totentries files; check what's going on")
                    exit()

                # os.chdir(jobname)

                #run as python nano_to_atto.py -l mu -y 2018 -m hadd -d singlemuonA -o /store/user/lpcrutgers/sthayil/pseudoaxions/atto_passTrigger -hd /store/user/lpcrutgers/sthayil/pseudoaxions/atto_hadded
                if os.path.isfile('haddfilelist'):
                    command='python scripts/haddnano.py /eos/uscms'+ops.haddDir+'/'+jobname+'.root '
                    for rootFile in rootFiles:
                        command+='/eos/uscms'+ops.outputDirectory+'/'+jobname+'/'+rootFile+' '
                    print(command)                    
                    os.system(command)
                    os.system('rm haddfilelist')

                # command='python ../scripts/haddnano.py '+jobname+'.root '
                # for rootFile in rootFiles: 
                #     os.system('xrdcp -s root://cmseos.fnal.gov/'+ops.outputDirectory+'/'+jobname+'/'+rootFile+' .')
                #     command+=rootFile+' '
                # os.system(command)
                # os.system('mv '+jobname+'.root ../'+ops.haddDir)
                # for rootFile in rootFiles: os.system('rm '+rootFile)

                tottotEntries=0
                print(totEntrFiles)
                for totEntrFile in totEntrFiles: 
                    os.system('xrdcp -s root://cmseos.fnal.gov/'+ops.outputDirectory+'/'+jobname+'/'+totEntrFile+' .')
                    with open(totEntrFile, 'r') as in_file: 
                        dataline=in_file.readline()
                        tottotEntries+=int(dataline)
                    os.system('rm '+totEntrFile)

                with open(jobname+'_totentries.txt', 'w') as entriesfile:
                    entriesfile.write(str(tottotEntries))
                #os.system('xrdcp '+jobname+'_totentries.txt root://cmseos.fnal.gov/'+ops.haddDir+'/'+jobname+'_totentries.txt')
                #print('xrdcp '+jobname+'_totentries.txt root://cmseos.fnal.gov/'+ops.haddDir+'/'+jobname+'_totentries.txt')
                #os.system('rm '+jobname+'_totentries.txt')

                # with open('../'+ops.haddDir+'/'+jobname+'_totentries.txt', 'w') as entriesfile:
                #     entriesfile.write(str(tottotEntries))
                # os.system('rm haddfilelist')

        else:
            print ("ERROR: Output directory is not an EOS path starting in /store/user/")
            exit()
    elif "hexcms" in hostname: 
        if not os.path.isdir(ops.haddDir): 
            print ("Making job directory "+ops.haddDir)
            os.system('mkdir '+ops.haddDir)

        if os.path.isdir(ops.outputDirectory+'/'+jobname):
            os.system('ls '+ops.outputDirectory+'/'+jobname+' > haddfilelist')
            with open('haddfilelist', 'r+') as in_file:
                rootFiles=[]
                totEntrFiles=[]
                for count, line in enumerate(in_file):
                    if ".root" in line: rootFiles.append(line.rstrip('\n'))
                    elif ".txt" in line: totEntrFiles.append(line.rstrip('\n'))
                if len(rootFiles)!=len(totEntrFiles): 
                    print ("#root files != #totentries files; check what's going on")
                    exit()

                os.system('rm haddfilelist')
                os.chdir(jobname)

                command='python ../scripts/haddnano.py '+jobname+'.root '
                for rootFile in rootFiles: 
                    #os.system('xrdcp -s root://cmseos.fnal.gov/'+ops.outputDirectory+'/'+jobname+'/'+rootFile+' .')
                    command+=ops.outputDirectory+'/'+jobname+'/'+rootFile+' '
                os.system(command)
                os.system('mv '+jobname+'.root ../'+ops.haddDir)
                #for rootFile in rootFiles: os.system('rm '+rootFile)
                tottotEntries=0
                for totEntrFile in totEntrFiles: 
                    #os.system('xrdcp -s root://cmseos.fnal.gov/'+ops.outputDirectory+'/'+jobname+'/'+totEntrFile+' .')
                    with open(ops.outputDirectory+'/'+jobname+'/'+totEntrFile, 'r') as in_file: 
                        dataline=in_file.readline()
                        tottotEntries+=int(dataline)
                    #os.system('rm '+totEntrFile)
                print ("\n\n"+str(tottotEntries))
                with open('../'+ops.haddDir+'/'+jobname+'_totentries.txt', 'w') as entriesfile:
                    entriesfile.write(str(tottotEntries))
        else:
            print ("ERROR: Output directory: "+ops.outputDirectory+'/'+jobname+" is not a valid path")
    exit()

#directory with log files--------------------------------------------------------------------
jobname = ops.year+'_'+ops.dataset+'_'+ops.lepton
if not os.path.isdir(jobname): 
    print ("Making job directory "+jobname)
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
        print ("ERROR: for output directory, specify an EOS path starting in /store/user/")
        exit()
    prefix="root://cmseos.fnal.gov/"
elif "hexcms" in hostname: 
    if os.path.isdir(ops.outputDirectory):
        if not ( (ops.outputDirectory).startswith("/") ): prefix=os.getcwd()+"/" #for hex, don't care about eos area, can save ouput locally
    else:
        print( "provided output directory does not exist, trying to create...")
        if not ( (ops.outputDirectory).startswith("/") ): prefix=os.getcwd()+"/"
        os.system('mkdir -p '+ops.outputDirectory)
        print( "created directory: "+prefix+ops.outputDirectory)
    if os.path.isdir(ops.outputDirectory+'/'+jobname): 
        print( "WARNING: directory named "+ops.outputDirectory+'/'+jobname+" already exists, files will be overwritten")
    else: os.system('mkdir -p '+ops.outputDirectory+'/'+jobname)
outputDir=prefix+ops.outputDirectory+'/'+jobname

#check accessibility of filelists-----------------------------------------------------------
if ".fnal.gov" in hostname:
    filelist = 'CMSSW_14_0_0/src/PhysicsTools/NanoAODTools/filelists/'+ops.year+'_'+ops.dataset+'.txt' 
    if not os.path.exists(filelist): 
        print( "Check that your filelist exists (format: filelists/year_dataset.txt)")
        exit
    f = open(filelist, "r")
    flines = f.readlines()

    if not (flines[0].startswith('/store/user/')): print( "Filelist items don't start in /store/user/; first item is: "+flines[0])
    f.close()
    
    #print #files/job, cd to jobdir, calculate num jobs------------------------------------------
    print( "Each job will run over up to "+str(ops.numFiles)+" NanoAODs")

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
    print ("Submitting jobs with condor_submit condorsubmit_nanotoatto_"+mytimestr+".jdl ...")
    #if ops.schedd is not None: os.system('condor_submit -name lpcschedd'+ops.schedd+'.fnal.gov condorsubmit_lhetominiaod_'+mytimestr+'.jdl')
    os.system('condor_submit condorsubmit_lhetominiaod_'+mytimestr+'.jdl')
    os.system('mv condorsubmit_lhetominiaod_'+mytimestr+'.jdl jdl_files/')

elif "hexcms" in hostname:
    directory = "CMSSW_14_0_0/src/PhysicsTools/NanoAODTools/filelists/"
    pattern = re.compile(rf"{ops.year}_{ops.dataset}(_\d+)?\.txt")
    matching_files = []

    print("Will run over these filelists:")
    for file_name in os.listdir(directory):
        if pattern.match(file_name):
            print(file_name)
            matching_files.append(file_name)

    numfiles = len(matching_files)

    maxfilesperfilelist=0
    numjobs=0
    for filelist in matching_files:
        if not os.path.exists(directory+filelist):
            print( "Check that your filelist exists (format: filelists/year_dataset.txt); provided: ",directory+filelist)
            exit
        f = open(directory+filelist, "r")
        flines = f.readlines()

        if not (flines[0].startswith('/cms/')): print( "Filelist items don't start in /cms/; first item is: "+flines[0])

        if len(flines)>maxfilesperfilelist: maxfilesperfilelist=len(flines)
        numjobsthislist = len(flines)//ops.numFiles + (len(flines)%ops.numFiles >0)
        print( "This filelist has "+str(len(flines))+" NanoAODs and will be split into "+str(numjobsthislist)+" jobs")
        numjobs+=numjobsthislist

    maxnumjobsperfilelist=maxfilesperfilelist//ops.numFiles + (maxfilesperfilelist%ops.numFiles >0)     
    os.chdir(jobname)
    os.system('cp /tmp/x509up_u775 ./x509up')
    if not os.path.isdir('logs'): os.system('mkdir logs')
    if not os.path.isdir('jdl_files'): os.system('mkdir jdl_files')

    #modify jdl---------------------------------------------------------------------------------
    mytime=datetime.now()
    mytimestr=mytime.strftime("%y%m%d%H")
    with open('../hex_condorsubmit_nanotoatto.jdl', 'r') as file :
        filedata = file.read()
    filedata = filedata.replace('00000000', ops.lepton) #ops.lepton ops.dataset ops.year outputDir ops.numFiles numjobs
    filedata = filedata.replace('11111111', ops.dataset)
    filedata = filedata.replace('22222222', ops.year)
    filedata = filedata.replace('33333333', outputDir)
    filedata = filedata.replace('44444444', str(ops.numFiles))
    filedata = filedata.replace('66666666', str(maxnumjobsperfilelist))
    filedata = filedata.replace('55555555', str(numjobs))
    with open('hex_condorsubmit_lhetominiaod_'+mytimestr+'.jdl', 'w') as file:
        file.write(filedata)
    
    #submit jobs-------------------------------------------------------------------------------
    os.system('cp ../hex_condorsubmit_nanotoatto.sh .')
    print ("Submitting jobs with condor_submit hex_condorsubmit_nanotoatto_"+mytimestr+".jdl ...")
    #if ops.schedd is not None: os.system('condor_submit -name lpcschedd'+ops.schedd+'.fnal.gov condorsubmit_lhetominiaod_'+mytimestr+'.jdl')
    os.system('condor_submit hex_condorsubmit_lhetominiaod_'+mytimestr+'.jdl')
    os.system('mv hex_condorsubmit_lhetominiaod_'+mytimestr+'.jdl jdl_files/')

    
#totalEvents mode-------------------------------------------------------------------------
