import sys, os, socket

#usage: python makeFilelist.py /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/ttPhiPS_M-250/ eosdirtocheck nameofexpectedfilelist
hostname = socket.gethostname()
filepath = sys.argv[1] 
if ".fnal.gov" in hostname: 
    os.system("eos root://cmseos.fnal.gov ls "+filepath+" > filelist")
    if "/eos/uscms/" in filepath: filepath = filepath.replace("/eos/uscms","")
elif "hexcms" in hostname: os.system("ls "+filepath+" > filelist")

with open(sys.argv[2], 'w') as out_file:
    with open('filelist', 'r') as in_file:
        for line in in_file:
            if "reports" in line: continue
            out_file.write(filepath + line.rstrip('\n') + '\n')

os.system("rm filelist")
