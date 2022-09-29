import os, sys

badfiles=[]
command='python ../scripts/haddnano.py '+sys.argv[4]+'_'+sys.argv[1]

for i in range(int(sys.argv[2]), int(sys.argv[3])):
    if i not in badfiles: command+=' NANOAOD_TwoProng_'+str(i)+'_Skim.root'

os.system(command)
