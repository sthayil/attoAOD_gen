import os, sys
import glob, argparse, socket, time
options = argparse.ArgumentParser(description="Find missing files")
options.add_argument("-i", "--inputDir", required=True, help="Dir to check")
options.add_argument("-o", "--outputFilename", required=True, help="name of output file (without .root)")
ops = options.parse_args()

filenums=[]
os.system("ls "+ops.inputDir+" > filelist")
dirfiles=open('filelist','r')
for line in dirfiles.readlines():
    if line.startswith('NANOAOD'): 
        filenum=int(filter(str.isdigit, line))
        filenums.append(filenum)
dirfiles.close()

print "num of files: ", len(filenums)
print "max filenum: ", max(filenums)

for j in range(0,int(max(filenums)/1000)+1):
    command='python scripts/haddnano.py '+str(j)+'_'+ops.outputFilename+'.root '
    for k in range(0,1000):
        if (1000*j)+k in filenums: command+=ops.inputDir+'/NANOAOD_TwoProng_'+str((1000*j)+k)+'_Skim.root '
    os.system(command)

command='python scripts/haddnano.py '+ops.outputFilename+'.root '
for j in range(0,int(max(filenums)/1000)+1): command+=str(j)+'_'+ops.outputFilename+'.root '
os.system(command)

for j in range(0,int(max(filenums)/1000)+1): os.system('rm '+str(j)+'_'+ops.outputFilename+'.root ')
os.system("rm filelist")
