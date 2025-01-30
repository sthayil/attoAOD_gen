from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import os, glob, sys, argparse, socket
from datetime import datetime
import json
from correctionlib import CorrectionSet

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
#240216: accomodate 2017 singlemuonB spl trigger reqs + add highptid req to goodmuon
#240918: copy 2016 triggers for 2016apv (not sure if legit); save 20160 into year branch for 2016APV
#241011: first pass at muon reco/id/iso/trig sfs; more stringent selection

with open('ScaleFactors_Muon_highPt_RECO_2018_schemaV2.json', 'r') as reco:
    correction_data_mu_reco = json.load(reco)
with open('ScaleFactors_Muon_highPt_IDISO_2018_schemaV2.json', 'r') as idiso:
    correction_data_mu_idiso = json.load(idiso)
with open('ScaleFactors_Muon_highPt_HLT_2018_schemaV2.json', 'r') as hlt:
    correction_data_mu_trig = json.load(hlt)

correction_set_mu_reco = CorrectionSet.from_json(correction_data_mu_reco)
correction_set_mu_idiso = CorrectionSet.from_json(correction_data_mu_idiso)
correction_set_mu_trig = CorrectionSet.from_json(correction_data_mu_trig)

# def get_mu_reco_scale_factor(abseta, pt, scale_factor_type):
#     correction = correction_set_reco['NUM_GlobalMuons_DEN_TrackerMuonProbes']
#     return correction.evaluate(abseta, pt, scale_factor_type)
# def get_mu_id_scale_factor(abseta, pt, scale_factor_type):
#     correction = correction_set_idiso['NUM_HighPtID_DEN_GlobalMuonProbes'] 
#     return correction.evaluate(abseta, pt, scale_factor_type)
# def get_mu_iso_scale_factor(abseta, pt, scale_factor_type):
#     correction = correction_set_idiso['NUM_probe_LooseRelTkIso_DEN_HighPtProbes'] 
#     return correction.evaluate(abseta, pt, scale_factor_type)
# def get_mu_hlt_scale_factor(abseta, pt, scale_factor_type):
#     correction = correction_set_trig['NUM_HLT_DEN_HighPtLooseRelIsoProbes'] 
#     return correction.evaluate(abseta, pt, scale_factor_type)

def get_mu_combined_scale_factor(abseta, pt, scale_factor_type):
    correction_reco = correction_set_reco['NUM_GlobalMuons_DEN_TrackerMuonProbes']
    correction_id   = correction_set_idiso['NUM_HighPtID_DEN_GlobalMuonProbes']
    correction_iso  = correction_set_idiso['NUM_probe_LooseRelTkIso_DEN_HighPtProbes'] 
    correction_hlt  = correction_set_trig['NUM_HLT_DEN_HighPtLooseRelIsoProbes'] 
    return correction_reco.evaluate(abseta, pt, scale_factor_type)*correction_id.evaluate(abseta, pt, scale_factor_type)*correction_iso.evaluate(abseta, pt, scale_factor_type)*correction_hlt.evaluate(abseta, pt, scale_factor_type)

def check_goodtwoprong(twoprongs):
    if not any(twoprong.pt > 20 and abs(twoprong.eta) < 1.3 for twoprong in twoprongs):
        return False

def check_goodmuon(muons):
    if not any(muon.pt > 52 and abs(muon.eta) < 2.4 and muon.highPtId == 2 for muon in muons):
        return False

def check_goodmuon(muons):
    good_muons = [muon for muon in muons if muon.pt > 52 and abs(muon.eta) < 2.4 and muon.highPtId == 2]
    if not good_muons:
        return False, None
    
    leading_muons = [muon for muon in good_muons if muon.tkRelIso < 0.1]
    if not leading_muons:
        return True, None
    
    leading_muon = max(leading_muons, key=lambda muon: muon.pt)
    return True, leading_muon

def check_goodjets(jets):
    if len(jets) < 3:
        return False
    goodjets = sum(1 for jet in jets if jet.pt > 30 and abs(jet.eta) < 2.5)
    return goodjets >= 3

selections=True #req twoprong, basic lep, pass trig

class attoAOD_ttw_mu(Module):
    def __init__(self, year="2018", mctype="0", attoVersion="241011"): 
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
        #self.out.branch("Muon_sf_nominal","F",lenVar="nMuon")
        self.out.branch("Muon_sf_nominal","F")
        self.out.branch("Muon_sf_systup", "F")
        self.out.branch("Muon_sf_systdown","F")
        self.out.branch("Muon_sf_stat","F")
        self.out.branch("Muon_sf_syst","F")
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
        if (self.year=="2018" or self.year=="2017"):
            if hasattr(event, "HLT_OldMu100"):
                if (trigger.Mu50 or trigger.OldMu100 or trigger.TkMu100): 
                    self.out.fillBranch("passTrigger", True)
            else:
                if (trigger.Mu50):
                    self.out.fillBranch("passTrigger", True)
        elif (self.year=="2016" or self.year=="2016APV") and (trigger.Mu50 or trigger.TkMu50): 
            self.out.fillBranch("passTrigger", True)
        else: 
            self.out.fillBranch("passTrigger", False)
            if selections: return False

        #event selection
        good_muon_exists, leading_iso_muon = check_goodmuon(muons)
        if selections:
            if not ( check_goodtwoprong(twoprongs) and check_goodmuon(muons) and check_goodjets(jets) ): return False

        #fill branches
        if self.year != "2016APV": self.out.fillBranch("year", int(self.year))
        else: self.out.fillBranch("year", int(20160))
        self.out.fillBranch("mcType", int(self.mctype))
        self.out.fillBranch("attoVersion", int(self.attoVersion))
        if leading_iso_muon:
            self.out.fillBranch("Muon_sf_nominal", get_mu_combined_scale_factor(abs(leading_iso_muon.eta), leading_iso_muon.pt, "nominal"))
            self.out.fillBranch("Muon_sf_systup",  get_mu_combined_scale_factor(abs(leading_iso_muon.eta), leading_iso_muon.pt, "systup"))
            self.out.fillBranch("Muon_sf_systdown",get_mu_combined_scale_factor(abs(leading_iso_muon.eta), leading_iso_muon.pt, "systdown"))
            self.out.fillBranch("Muon_sf_stat",    get_mu_combined_scale_factor(abs(leading_iso_muon.eta), leading_iso_muon.pt, "stat"))
            self.out.fillBranch("Muon_sf_syst",    get_mu_combined_scale_factor(abs(leading_iso_muon.eta), leading_iso_muon.pt, "syst"))
        else:
            self.out.fillBranch("Muon_sf_nominal",-99.0)
            self.out.fillBranch("Muon_sf_systup",-99.0)
            self.out.fillBranch("Muon_sf_systdown",-99.0)
            self.out.fillBranch("Muon_sf_stat",-99.0)
            self.out.fillBranch("Muon_sf_syst",-99.0)

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
mu_2018_sig_M3000  =  lambda: attoAOD_ttw_mu(mctype="17")
mu_2018_sig_M3000_tcoupling1_Wb1Wb1  =  lambda: attoAOD_ttw_mu(mctype="18")
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
mu_2016APV_singlemuonB = lambda: attoAOD_ttw_mu(year="2016APV")
mu_2016APV_singlemuonC = lambda: attoAOD_ttw_mu(year="2016APV")
mu_2016APV_singlemuonD = lambda: attoAOD_ttw_mu(year="2016APV")
mu_2016APV_singlemuonE = lambda: attoAOD_ttw_mu(year="2016APV")
mu_2016APV_singlemuonF = lambda: attoAOD_ttw_mu(year="2016APV")
mu_2016APV_singlemuonF2= lambda: attoAOD_ttw_mu(year="2016APV")
mu_2016APV_singlemuonG2= lambda: attoAOD_ttw_mu(year="2016APV")
mu_2016APV_singlemuonH2= lambda: attoAOD_ttw_mu(year="2016APV")
mu_2016APV_ttjets =      lambda: attoAOD_ttw_mu(year="2016APV",mctype="20")
mu_2016APV_wjetstolnu =  lambda: attoAOD_ttw_mu(year="2016APV",mctype="21")
mu_2016APV_dyjetstoll =  lambda: attoAOD_ttw_mu(year="2016APV",mctype="22")
mu_2016APV_sig_M1000 =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="11")
mu_2016APV_sig_M500  =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="12")
mu_2016APV_sig_M250  =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="13")
mu_2016APV_sig_M4000 =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="14")
mu_2016APV_sig_M2000  =  lambda: attoAOD_ttw_mu(year="2016APV",mctype="15")
mu_2016APV_sig_M750  =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="16")
