from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import os, glob, sys, argparse, socket
from datetime import datetime

#0: /SingleMuon/Run2018B-UL2018_MiniAODv2-v2/MINIAOD
#11:/eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-1000/
#12:/eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-500/
#13:/eos/uscms/store/user/lpcrutgers/sthayil/pseudoaxions/nano/ttPhiPS_M-250/
#21: /WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM  
#    ( NO LONGER /TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM )
#22: /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM

#220929: added attoversion branch, changed mctype def to add wjets
#221007: require that events pass trigger to be saved into atto (to reduce filesizes)

class attoAOD_ttw(Module):
    def __init__(self, finalStateLepton="mu", year="2018", mctype="0", attoVersion="221007"): #el/.../0,21,22
        self.finalStateLepton = finalStateLepton
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
        self.out.branch("Muon_scaleFactor","F",lenVar="nMuon")
        self.out.branch("Electron_scaleFactor","F",lenVar="nMuon")
        self.out.branch("year","I")
        self.out.branch("mcType","I")
        self.out.branch("passTrigger","O")
        self.out.branch("attoVersion","I")
        #self.out.branch("evProcessed","I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        flags = Object(event, "Flag")
        trigger = Object(event, "HLT")
        jets = Collection(event, "Jet")
        twoprongs = Collection(event, "TwoProng")

        if self.finalStateLepton == 'mu' :  leptons = Collection(event, "Muon")
        elif self.finalStateLepton == 'el' :  leptons = Collection(event, "Electron")
        else: return False

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
            if not (pass_filters): return False

        #trigger
        if self.finalStateLepton == 'mu' and (trigger.Mu50 or trigger.OldMu100 or trigger.TkMu100): self.out.fillBranch("passTrigger", True)
        elif self.finalStateLepton == 'el' and (trigger.Ele50_CaloIdVT_GsfTrkIdT_PFJet165 or trigger.Photon200): self.out.fillBranch("passTrigger", True)
        else: 
            self.out.fillBranch("passTrigger", False)
            return False #temporary

        #event selection
        goodTwoprong=False
        if len(twoprongs)<1: return False
        for twoprong in twoprongs:
            if twoprong.pt>20 and abs(twoprong.eta)<2.5: goodTwoprong=True
        if goodTwoprong==False: return False

        goodLepton=False
        if len(leptons)<1: return False
        for lepton in leptons:
            if lepton.pt>52 and abs(lepton.eta)<2.5: 
                goodLepton = True
        if goodLepton==False: return False

        #fill branches
        self.out.fillBranch("year", int(self.year))
        self.out.fillBranch("mcType", int(self.mctype))
        self.out.fillBranch("attoVersion", int(self.attoVersion))
        #musf=[1 for x in range(len(muons))]
        #self.out.fillBranch("Muon_scaleFactor",musf)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
attoAOD_ttw_mu_2018_data =        lambda: attoAOD_ttw()
attoAOD_ttw_mu_2018_ttjets =       lambda: attoAOD_ttw(mctype="20")
attoAOD_ttw_mu_2018_wjetstolnu = lambda: attoAOD_ttw(mctype="21")
#attoAOD_ttw_mu_2018_tttosemilep = lambda: attoAOD_ttw(mctype="21")
attoAOD_ttw_mu_2018_dyjetstoll =  lambda: attoAOD_ttw(mctype="22")
attoAOD_ttw_el_2018_data =        lambda: attoAOD_ttw(finalStateLepton="el")
attoAOD_ttw_el_2018_ttjets =       lambda: attoAOD_ttw(finalStateLepton="el",mctype="20")
attoAOD_ttw_el_2018_wjetstolnu = lambda: attoAOD_ttw(finalStateLepton="el",mctype="21")
#attoAOD_ttw_el_2018_tttosemilep = lambda: attoAOD_ttw(finalStateLepton="el",mctype="21")
attoAOD_ttw_el_2018_dyjetstoll =  lambda: attoAOD_ttw(finalStateLepton="el",mctype="22")
