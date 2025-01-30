import sys, os, socket

#usage: python makeFilelist.py /store/user/lpcrutgers/sthayil/pseudoaxions/ttPhiPS_M-250/ eosdirtocheck nameofexpectedfilelist
#new usage: python3 makeFilelist.py /cms/twoprong/johnpaul/crab_Nov8/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/DYJets2017/241127_210125/ filelists/2017_dyjetstoll.txt
hostname = socket.gethostname()
filepath = sys.argv[1] 
if ".fnal.gov" in hostname: 
    os.system("eos root://cmseos.fnal.gov ls "+filepath+" > filelist")
    if "/eos/uscms/" in filepath: filepath = filepath.replace("/eos/uscms","")

    with open(sys.argv[2], 'w') as out_file:
        with open('filelist', 'r') as in_file:
            for line in in_file:
                if "reports" in line: continue
                out_file.write(filepath + line.rstrip('\n') + '\n')
        
    os.system("rm filelist")


elif "hexcms" in hostname:
    root_files = [f for f in os.listdir(filepath) if f.endswith(".root")]
    directories = sorted([d for d in os.listdir(filepath) if os.path.isdir(os.path.join(filepath, d))])

    if not root_files and directories:
        output_basename = sys.argv[2].rsplit('.txt', 1)[0]
        
        print("No .root files found. The following directories exist:")
        for dircnt, directory in enumerate(directories):
            print(directory)

            full_directory_path = os.path.join(filepath, directory)
            os.system(f"ls {full_directory_path} > filelist")
            output_filename = f"{output_basename}_{dircnt}.txt"

            with open(output_filename, 'w') as out_file:
                with open('filelist', 'r') as in_file:
                    for line in in_file:
                        if "reports" in line or "log" in line:
                            continue
                        out_file.write(filepath + "/" + directory + "/" + line.rstrip('\n') + '\n')

            os.system("rm filelist")
    
    elif root_files:
        os.system("ls "+filepath+" > filelist")

        output_basename = sys.argv[2].rsplit('.txt', 1)[0] #hack
        output_filename = f"{output_basename}_0.txt" #hack
        with open(output_filename, 'w') as out_file:
            with open('filelist', 'r') as in_file:
                for line in in_file:
                    if "reports" in line: continue
                    out_file.write(filepath + line.rstrip('\n') + '\n')
        
        os.system("rm filelist")
