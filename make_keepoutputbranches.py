import os

#usage: python make_keepoutputbranches.py

with open('keepoutputbranches_mu.txt', 'w') as out_file_1:
    with open('keepoutputbranches.txt', 'r') as in_file:
        for line in in_file:
            if line.strip() == '': continue
            out_file_1.write(line.rstrip('\n') + '\n')
    out_file_1.write('keep Muon*\n')
    out_file_1.write('keep HLT_Mu50\n')
    out_file_1.write('keep HLT_OldMu100\n')
    out_file_1.write('keep HLT_TkMu100\n')

with open('keepoutputbranches_el.txt', 'w') as out_file_2:
    with open('keepoutputbranches.txt', 'r') as in_file:
        for line in in_file:
            if line.strip() == '': continue
            out_file_2.write(line.rstrip('\n') + '\n')
    out_file_2.write('keep Electron*\n')
    out_file_2.write('keep HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165\n')
    out_file_2.write('keep HLT_Photon200\n')


