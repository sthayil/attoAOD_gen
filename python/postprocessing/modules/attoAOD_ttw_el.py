from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import os, glob, sys, argparse, socket
from datetime import datetime

#0:  //EGamma/Run2018*-UL2018_MiniAODv2-v*/MINIAOD
#11: /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-1000/
#12: /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-500/
#13: /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-250/
#14: /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-4000/2023-04-24-11-00-18/fv1p5-9-5ce1_bv1p2-0-1531/
#15: /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-2000/2023-04-24-10-59-20/fv1p5-9-5ce1_bv1p2-0-1531/
#16: /eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-750/2023-04-24-11-00-43/fv1p5-9-5ce1_bv1p2-0-1531/
#20: /TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM
#21: /WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM  
#22: /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM

#220929: added attoversion branch, changed mctype def to add wjets
#221007: require that events pass trigger to be saved into atto (to reduce filesizes)
#221013: split el and mu modules
#230412: include signal MCs in options [with all selections (twoprong, lep and trig)]
#230529: introduce 'selections' boolean

selections=False

class attoAOD_ttw_el(Module):
    def __init__(self, year="2018", mctype="0", attoVersion="230529"): 
        self.year = year
        self.mctype = mctype
        self.attoVersion = attoVersion
        pass

    def mygetattr(self, my_obj, my_branch, default_bool):
        try: getattr(my_obj, my_branch)
        except RuntimeError: return default_bool
        else: return getattr(my_obj, my_branch)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        #self.out.branch("Electron_scaleFactor","F",lenVar="nElectron")
        self.out.branch("year","I")
        self.out.branch("mcType","I")
        self.out.branch("passTrigger","O")
        self.out.branch("attoVersion","I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        flags = Object(event, "Flag")
        trigger = Object(event, "HLT")
        jets = Collection(event, "Jet")
        twoprongs = Collection(event, "TwoProng")
        electrons = Collection(event, "Electron")

        # baseline filtering
        if self.mctype=="0":
            pass_filters = (
                self.mygetattr(flags, 'goodVertices', True)
                and self.mygetattr(flags, 'HBHENoiseFilter', True)
                and self.mygetattr(flags, 'HBHENoiseIsoFilter', True)
                and self.mygetattr(flags, 'EcalDeadCellTriggerPrimitiveFilter', True)
                and self.mygetattr(flags, 'BadPFMuonFilter', True)
                and self.mygetattr(flags, 'BadChargedCandidateFilter', True)
                and self.mygetattr(flags, 'ecalBadCalibFilter', True)
                and self.mygetattr(flags, 'globalSuperTightHalo2016Filter', True)
                and self.mygetattr(flags, 'eeBadScFilter', True)
            )
            if not (pass_filters): 
                return False

        #trigger
        if (trigger.Ele50_CaloIdVT_GsfTrkIdT_PFJet165 or trigger.Photon200): 
            self.out.fillBranch("passTrigger", True)
        else: 
            self.out.fillBranch("passTrigger", False)
            if selections: return False 

        #event selection
        if selections:
            goodTwoprong=False
            if len(twoprongs)<1: return False
            for twoprong in twoprongs:
                if twoprong.pt>20 and abs(twoprong.eta)<2.5: 
                    goodTwoprong=True
                    break
            if goodTwoprong==False: return False

            goodElectron=False
            if len(electrons)<1: return False
            for electron in electrons:
                if electron.pt>52 and abs(electron.eta)<2.5: 
                    goodElectron = True
                    break
            if goodElectron==False: return False

        #fill branches
        self.out.fillBranch("year", int(self.year))
        self.out.fillBranch("mcType", int(self.mctype))
        self.out.fillBranch("attoVersion", int(self.attoVersion))
        #elsf=[1 for x in range(len(leptons))]
        #self.out.fillBranch("Electron_scaleFactor",elsf)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
el_2018_data =        lambda: attoAOD_ttw_el()
el_2018_egammaA =     lambda: attoAOD_ttw_el()
el_2018_egammaB =     lambda: attoAOD_ttw_el()
el_2018_egammaC =     lambda: attoAOD_ttw_el()
el_2018_egammaD =     lambda: attoAOD_ttw_el()
el_2018_ttjets =      lambda: attoAOD_ttw_el(mctype="20")
el_2018_wjetstolnu =  lambda: attoAOD_ttw_el(mctype="21")
el_2018_dyjetstoll =  lambda: attoAOD_ttw_el(mctype="22")
el_2018_sig_M1000 =  lambda: attoAOD_ttw_el(mctype="11")
el_2018_sig_M500  =  lambda: attoAOD_ttw_el(mctype="12")
el_2018_sig_M250  =  lambda: attoAOD_ttw_el(mctype="13")
el_2018_sig_M4000 =  lambda: attoAOD_ttw_el(mctype="14")
el_2018_sig_M2000  = lambda: attoAOD_ttw_el(mctype="15")
el_2018_sig_M750  =  lambda: attoAOD_ttw_el(mctype="16")
