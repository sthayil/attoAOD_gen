from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import os, glob, sys, argparse, socket
from datetime import datetime

#0:  /SingleMuon/Run2018*-UL2018_MiniAODv2-v*/MINIAOD
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
#221209: undo trig req
#230227: undo all selections (no twoprong, lep or trig)
#230410: redo all selections (twoprong, lep and trig)
#230412: include signal MCs in options [with all selections (twoprong, lep and trig)]
#230529: introduce 'selections' boolean
#240215: change triggers to accomodate 2017, 2016 mu

selections=True #req twoprong, basic lep, pass trig

class attoAOD_ttw_mu(Module):
    def __init__(self, year="2018", mctype="0", attoVersion="240215"): 
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
        #self.out.branch("Muon_scaleFactor","F",lenVar="nMuon")
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
        muons = Collection(event, "Muon")

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
        if (self.year=="2018" or self.year=="2017") and (trigger.Mu50 or trigger.OldMu100 or trigger.TkMu100): 
            self.out.fillBranch("passTrigger", True)
        elif self.year=="2016" and (trigger.Mu50 or trigger.TkMu50): 
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

            goodMuon=False
            if len(muons)<1: return False
            for muon in muons:
                if muon.pt>52 and abs(muon.eta)<2.4: 
                    goodMuon = True
                    break
            if goodMuon==False: return False

        #fill branches
        self.out.fillBranch("year", int(self.year))
        self.out.fillBranch("mcType", int(self.mctype))
        self.out.fillBranch("attoVersion", int(self.attoVersion))
        #musf=[1 for x in range(len(leptons))]
        #self.out.fillBranch("Muon_scaleFactor",musf)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
mu_2018_singlemuonA = lambda: attoAOD_ttw_mu()
mu_2018_singlemuonB = lambda: attoAOD_ttw_mu()
mu_2018_singlemuonC = lambda: attoAOD_ttw_mu()
mu_2018_singlemuonD = lambda: attoAOD_ttw_mu()
mu_2018_ttjets =      lambda: attoAOD_ttw_mu(mctype="20")
mu_2018_wjetstolnu =  lambda: attoAOD_ttw_mu(mctype="21")
mu_2018_dyjetstoll =  lambda: attoAOD_ttw_mu(mctype="22")
mu_2018_sig_M1000 =   lambda: attoAOD_ttw_mu(mctype="11")
mu_2018_sig_M500  =   lambda: attoAOD_ttw_mu(mctype="12")
mu_2018_sig_M250  =   lambda: attoAOD_ttw_mu(mctype="13")
mu_2018_sig_M4000 =   lambda: attoAOD_ttw_mu(mctype="14")
mu_2018_sig_M2000   = lambda: attoAOD_ttw_mu(mctype="15")
mu_2018_sig_M750  =   lambda: attoAOD_ttw_mu(mctype="16")
mu_2017_singlemuonB = lambda: attoAOD_ttw_mu(year="2017")
mu_2017_singlemuonC = lambda: attoAOD_ttw_mu(year="2017")
mu_2017_singlemuonD = lambda: attoAOD_ttw_mu(year="2017")
mu_2017_singlemuonE = lambda: attoAOD_ttw_mu(year="2017")
mu_2017_singlemuonF = lambda: attoAOD_ttw_mu(year="2017")
mu_2017_ttjets =      lambda: attoAOD_ttw_mu(year="2017",mctype="20")
mu_2017_wjetstolnu =  lambda: attoAOD_ttw_mu(year="2017",mctype="21")
mu_2017_dyjetstoll =  lambda: attoAOD_ttw_mu(year="2017",mctype="22")
mu_2017_sig_M1000 =   lambda: attoAOD_ttw_mu(year="2017",mctype="11")
mu_2017_sig_M500  =   lambda: attoAOD_ttw_mu(year="2017",mctype="12")
mu_2017_sig_M250  =   lambda: attoAOD_ttw_mu(year="2017",mctype="13")
mu_2017_sig_M4000 =   lambda: attoAOD_ttw_mu(year="2017",mctype="14")
mu_2017_sig_M2000  =  lambda: attoAOD_ttw_mu(year="2017",mctype="15")
mu_2017_sig_M750  =   lambda: attoAOD_ttw_mu(year="2017",mctype="16")
mu_2016_singlemuonB = lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonC = lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonD = lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonE = lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonF = lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonF2= lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonG2= lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonH2= lambda: attoAOD_ttw_mu(year="2016")
mu_2016_ttjets =      lambda: attoAOD_ttw_mu(year="2016",mctype="20")
mu_2016_wjetstolnu =  lambda: attoAOD_ttw_mu(year="2016",mctype="21")
mu_2016_dyjetstoll =  lambda: attoAOD_ttw_mu(year="2016",mctype="22")
mu_2016_sig_M1000 =   lambda: attoAOD_ttw_mu(year="2016",mctype="11")
mu_2016_sig_M500  =   lambda: attoAOD_ttw_mu(year="2016",mctype="12")
mu_2016_sig_M250  =   lambda: attoAOD_ttw_mu(year="2016",mctype="13")
mu_2016_sig_M4000 =   lambda: attoAOD_ttw_mu(year="2016",mctype="14")
mu_2016_sig_M2000  =  lambda: attoAOD_ttw_mu(year="2016",mctype="15")
mu_2016_sig_M750  =   lambda: attoAOD_ttw_mu(year="2016",mctype="16")
