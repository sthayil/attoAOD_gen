#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh 

pwd
date

printf "\n\n"
printf "Parameters passed:------------------------------------------------------------------------------------------"
printf "\n$1 : batch number"
printf "\n$2 : lepton (el/mu)"
printf "\n$3 : dataset"
printf "\n$4 : year (2016/2016APV/2017/2018)"
printf "\n$5 : outputDir to xrdcp to (EOS area for cmslpc)"
printf "\n$6 : (max) number of files to run over in this job"
printf "\n$7 : max number of jobs per filelist \n\n"
thisFilelist=$(( $1 / $7 ))
thisBatch=$(( $1 % $7 ))
printf "\n$thisFilelist : will run over filelist filelists/$4_$3_$thisFilelist.txt"
printf "\n$thisBatch    : will use batch number $thisBatch to run over filelist lines $thisBatch*$6 to ($thisBatch+1)*$6 \n\n"

mv  hex_CMSSW_14_0_0.tar.gz  CMSSW_14_0_0.tar.gz
tar -xf CMSSW_14_0_0.tar.gz
rm CMSSW_14_0_0.tar.gz
cd CMSSW_14_0_0/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`
echo $CMSSW_BASE

cd PhysicsTools/NanoAODTools
#cp ${_CONDOR_SCRATCH_DIR}/$4_$3.txt .
#cp ${_CONDOR_SCRATCH_DIR}/keepoutputbranches_$2.txt .

printf "\n\n"
printf "Current dir:------------------------------------------------------------------------------------------------\n"
pwd
date
ls -alh

python3 scripts/create_attoAOD.py -m $6 --batch $thisBatch -f filelists/$4_$3_$thisFilelist.txt -l $2 -I PhysicsTools.NanoAODTools.postprocessing.attoAOD_ttw_$2 $2_$4_$3 --bo keepoutputbranches_$2.txt


date
mv "$thisBatch"_$4_$3.root $1_$4_$3.root
cp totentries.txt $1_totentries.txt
printf "Current contents:------------------------------------------------------------------------------------------------\n"
ls -alhtr
xrdcp $1_$4_$3.root $5
xrdcp $1_totentries.txt $5
#rm $1_$2_$3_$4.root
