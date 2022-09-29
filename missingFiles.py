import sys, os

#usage: python missingfiles.py /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/ttPhiPS_M-250/eosdirtocheck NumOfExpectedFilesEg10

#os.system("eos root://cmseos.fnal.gov ls "+sys.argv[1]+" > filelist")
os.system("ls "+sys.argv[1]+" > filelist")
filenums=[]
dirfiles=open('filelist','r')
for line in dirfiles.readlines():
    if line.startswith('NANOAOD'): #to check for miniAOD_*.root, for example
        filenum=int(filter(str.isdigit, line))
        filenums.append(filenum)
dirfiles.close()

print "num of files: ", len(filenums)

for i in range(0, int(sys.argv[2])):
    if i not in filenums: print str(i)+",",

os.system("rm filelist")
