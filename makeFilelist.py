import sys, os, socket

#usage: python makeFilelist.py /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/ttPhiPS_M-250/ eosdirtocheck nameofexpectedfilelist
hostname = socket.gethostname()
if ".fnal.gov" in hostname: os.system("eos root://cmseos.fnal.gov ls "+sys.argv[1]+" > filelist")
elif "hexcms" in hostname: os.system("ls "+sys.argv[1]+" > filelist")

with open(sys.argv[2], 'w') as out_file:
    with open('filelist', 'r') as in_file:
        for line in in_file:
            out_file.write(sys.argv[1] + line.rstrip('\n') + '\n')

os.system("rm filelist")
