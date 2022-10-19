#!/bin/bash
#! /bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a tcsh script, use .csh instead of .sh
export SCRAM_ARCH=slc7_amd64_gcc700

pwd

eval `scramv1 project CMSSW CMSSW_10_6_20`
cd CMSSW_10_6_20/src
eval `scramv1 runtime -sh`
git clone https://github.com/sthayil/attoAOD_gen.git PhysicsTools/NanoAODTools
scramv1 b
eval `scramv1 runtime -sh`
cd PhysicsTools/NanoAODTools
cp ${_CONDOR_SCRATCH_DIR}/$3_$4.txt .
cp ${_CONDOR_SCRATCH_DIR}/keepoutputbranches_$2.txt .

printf "\n\n"
printf "Current dir:------------------------------------------------------------------------------------------------"
pwd
ls

printf "\n\n"
printf "Parameters passed:------------------------------------------------------------------------------------------"
printf "\n$1 : batch number"
printf "\n$2 : lepton (el/mu)"
printf "\n$3 : dataset"
printf "\n$4 : year (2016/2017/2018)"
printf "\n$5 : outputDir to xrdcp to (EOS area for cmslpc)"
printf "\n$6 : number of files to run over in this job\n\n"

python scripts/create_attoAOD.py -m $6 --batch $1 -f $3_$4.txt -l $2 -I PhysicsTools.NanoAODTools.postprocessing.modules.attoAOD_ttw_$2 $2_$4_$3 --bo keepoutputbranches_$2.txt

cp totentries.txt $1_totentries.txt
xrdcp $1_$2_$3_$4.root $5
xrdcp $1_totentries.txt $5
#rm $1_$2_$3_$4.root
