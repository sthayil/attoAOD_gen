import sys, os

#usage: python missingfiles.py /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/ttPhiPS_M-250/eosdirtocheck NumOfExpectedFilesEg10

import glob, argparse, socket
options = argparse.ArgumentParser(description="Find missing files")
options.add_argument("-m", "--mode", required=True, help="eos or local?", choices=['eos','local'])
options.add_argument("-i", "--inputDir", required=True, help="Dir to check")
ops = options.parse_args()

filenums=[]

if ops.mode == 'eos': os.system("eos root://cmseos.fnal.gov ls "+ops.inputDir+" > filelist")
elif ops.mode == 'local': os.system("ls "+ops.inputDir+" > filelist")

dirfiles=open('filelist','r')
for line in dirfiles.readlines():
    if line.startswith('NANOAOD'): #to check for miniAOD_*.root, for example
        filenum=int(filter(str.isdigit, line))
        filenums.append(filenum)
dirfiles.close()

print "num of files: ", len(filenums)
print "max filenum: ", max(filenums)

for i in range(0, max(filenums)+1):
    if i not in filenums: print str(i)+",",

os.system("rm filelist")
