from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import os, glob, sys, argparse, socket
from datetime import datetime
import json
from correctionlib import CorrectionSet
from functools import reduce
import random

#24-11-16: these need to change
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
#241015: split up muon sfs, remove stat, syst branches
#241116: add in new mass points, change mctypes, implement sfs only for signal mc
#241119: major additions for 2018 btag sfs; fix bad branches for 2017, 2016 trig
#241202: add muon, btag sfs to 2018 standard model mc
#241214: fix muon reco sf to use p instead of pt
#241216: add separate b-tagging efficiencies file for etaprime2018 (uses only M1000). also, bad trigger boolean from 241119 is now fixed.
#250129: add muon reco/id/iso/hlt sfs and b-tagging sfs for 2017, 2016, 2016apv. AND turn selections ON.
#250203: for 2016apv, 2016, 2017 non-signal mc: add b-tagging efficiency files. for mu reco/id/ido/hlt, include correction files inside mu sf function. for muon pts, introduce factor of leading_iso_muon.tunepRelPt*
#250217: turn selections OFF
#250405: changed sfs files filepaths; add muon momentum resolution smearing and systematics; add muon rochester corrections
#250406: add pileup reweighting SFs
#250429: modify muon p reso systematics to operate on the smeared p
#250515: add jet sfs and corrs (jec, jer) 
#!!!!!!!!!!!NOT DONE also for b-tagging signal mc sfs, make efficiency plots for each mass point.
version="250515"
print("AttoAOD version: ",version)
selections=True #req twoprong, basic lep, pass trig

ROOT.gSystem.Load("MuonRocCor/RoccoR.so")  #pre-compiled library

def correct_jet(year,mctype,run,rawFactor,pt,eta,mass,area,rho,genjetpt):
    #pT of the jet before specific correction (for JER and uncertainties: after all corrections
    scaletag="L1L2L3Res"
    resotag="PtResolution"
    resosftag="ScaleFactor" #nom, up, down
    correl_uncert=  ["Absolute", "FlavorQCD", "BBEC1", "EC2", "HF", "RelativeBal", "Total"]
    uncorrel_uncert=["Absolute", "BBEC1", "EC2", "HF", "RelativeSample"]
    
    if year=="2016APV":
        year="2016preVFP"
        yearforuncert="2016"
        if mctype==0: jectag="Summer19UL16APV_Run"+run+"_V7_DATA"
        else:
            jectag="Summer19UL16APV_V7_MC"
            jertag="Summer20UL16APV_JRV3_MC"

    elif year=="2016":
        year="2016postVFP"
        yearforuncert="2016"
        if mctype==0: jectag="Summer19UL16_RunFGH_V7_DATA"
        else:
            jectag="Summer19UL16_V7_MC"
            jertag="Summer20UL16_JRV3_MC"

    elif year=="2017":
        yearforuncert="2017"
        if mctype==0: jectag="Summer19UL17_Run"+run+"_V5_DATA"
        else:
            jectag="Summer19UL17_V5_MC"
            jertag="Summer19UL17_JRV2_MC"

    elif year=="2018":
        yearforuncert="2018"
        if mctype==0: jectag="Summer19UL18_Run"+run+"_V5_DATA"
        else:
            jectag="Summer19UL18_V5_MC"
            jertag="Summer19UL18_JRV2_MC"

    correction_file = "POG/JME/"+year+"_UL/jet_jerc.json.gz"
    correction_set  = CorrectionSet.from_file(correction_file)
    
    pt_raw = (1 - rawFactor)*pt

    correl_uncert_vals=[]
    uncorrel_uncert_vals=[]
    
    #evaluate corrections and uncertainties
    correction = correction_set.compound[jectag+"_"+scaletag+"_AK4PFchs"]
    corrval = correction.evaluate(area, eta, pt_raw, rho)
    pt_corr = pt_raw*corrval

    if mctype!=0:
        correction = correction_set[jertag+"_"+resosftag+"_AK4PFchs"]
        reso_basesf = correction.evaluate(eta,"nom")
        reso_basesf_up = correction.evaluate(eta,"up")
        reso_basesf_down = correction.evaluate(eta,"down")

        correction = correction_set[jertag+"_"+resotag+"_AK4PFchs"]
        ptreso = correction.evaluate(eta,pt_corr,rho)

        reso_sf = calculate_new_jet_resolution_sf(reso_basesf,pt_corr,genjetpt,ptreso)
        reso_sf_up = calculate_new_jet_resolution_sf(reso_basesf_up,pt_corr,genjetpt,ptreso)
        reso_sf_down = calculate_new_jet_resolution_sf(reso_basesf_down,pt_corr,genjetpt,ptreso)
        
        for uncert in correl_uncert:
            correction = correction_set[jectag+"_Regrouped_"+uncert+"_AK4PFchs"]
            correl_uncert_vals.append( correction.evaluate(eta,pt_corr) )
            
        for uncert in uncorrel_uncert:
            correction = correction_set[jectag+"_Regrouped_"+uncert+"_"+yearforuncert+"_AK4PFchs"]
            uncorrel_uncert_vals.append( correction.evaluate(eta,pt_corr) )

        return(corrval,pt_corr, reso_sf,reso_sf_up,reso_sf_down, correl_uncert_vals,uncorrel_uncert_vals)

    else: #for data
        return(corrval,pt_corr)

def calculate_new_jet_resolution_sf(reso_basesf,pt_corr,genjetpt,ptreso):
    if genjetpt<0: #stochastic smearing
        if max(reso_basesf**2 - 1, 0) == 0: reso_sf=1
        else:
            gauss_random = random.gauss(0, ptreso)
            reso_sf = 1 + (gauss_random * math.sqrt(reso_basesf**2 - 1))
    else: #scaling
        reso_sf = 1 + ( ((reso_basesf-1)*(pt_corr-genjetpt)) / pt_corr)
        if reso_sf<0: reso_sf=0
    return(reso_sf)
    
def get_pileup_weight(year,nTrueInt):
    if year=="2016APV": year="2016preVFP"
    elif year=="2016":  year="2016postVFP"
    correction_file = "POG/LUM/"+year+"_UL/puWeights.json.gz"
    pu_weight_corr = list(CorrectionSet.from_file(correction_file).values())[0] #this works, but just in case, keys are: Collisions16_UltraLegacy_goldenJSON, Collisions17_UltraLegacy_goldenJSON, Collisions18_UltraLegacy_goldenJSON
    pu_weight = pu_weight_corr.evaluate( nTrueInt, 'nominal' )
    pu_weight_up = pu_weight_corr.evaluate( nTrueInt, 'up' )
    pu_weight_down = pu_weight_corr.evaluate( nTrueInt, 'down' )
    return(pu_weight,pu_weight_up,pu_weight_down)
    
def get_rochester_mu_corrections(year,mctype,muon,genPt):
    if year=="2016APV": roccor_file="MuonRocCor/RoccoR2016aUL.txt" 
    elif year=="2016":  roccor_file="MuonRocCor/RoccoR2016bUL.txt" 
    elif year=="2017":  roccor_file="MuonRocCor/RoccoR2017UL.txt" 
    elif year=="2018":  roccor_file="MuonRocCor/RoccoR2018UL.txt" 
    roccor = ROOT.RoccoR(roccor_file)

    if mctype==0: #data
        roc_sf  = roccor.kScaleDT(int(muon.charge), muon.pt*muon.tunepRelPt, muon.eta, muon.phi, 0, 0)
        roc_err = roccor.kScaleDTerror(int(muon.charge), muon.pt*muon.tunepRelPt, muon.eta, muon.phi)
    else:
        if genPt != -999: #matched mc muon
            roc_sf  = roccor.kSpreadMC(int(muon.charge), muon.pt*muon.tunepRelPt, muon.eta, muon.phi, genPt, 0, 0)
            roc_err = roccor.kSpreadMCerror(int(muon.charge), muon.pt*muon.tunepRelPt, muon.eta, muon.phi, genPt)
        else: #unmatched mc muon
            nrandom = random.random()
            roc_sf  = roccor.kSmearMC(int(muon.charge), muon.pt*muon.tunepRelPt, muon.eta, muon.phi, int(muon.nTrackerLayers), nrandom, 0, 0)
            roc_err = roccor.kSmearMCerror(int(muon.charge), muon.pt*muon.tunepRelPt, muon.eta, muon.phi, int(muon.nTrackerLayers), nrandom)
    return(roc_sf, roc_err)
            
def muon_momentum_reso_sigma(year, p, eta):
    if abs(eta) < 1.2:
        if year=="2016APV": return 0.0112 + 6.87e-05 * p - 3.88e-08 * p**2 + 9.03e-12 * p**3
        elif year=="2016":  return 0.0102 + 6.77e-05 * p - 3.72e-08 * p**2 + 8.53e-12 * p**3
        elif year=="2017":  return 0.0104 + 6.11e-05 * p - 3.31e-08 * p**2 + 6.73e-12 * p**3
        elif year=="2018":  return 0.0108 + 5.93e-05 * p - 3.08e-08 * p**2 + 6.04e-12 * p**3
    elif 1.2 <= abs(eta) < 2.4:   
        if year=="2016APV": return 0.013  + 6.93e-05 * p - 3.46e-08 * p**2 + 7.72e-12 * p**3
        elif year=="2016":  return 0.0129 + 6.48e-05 * p - 3.04e-08 * p**2 + 6.63e-12 * p**3
        elif year=="2017":  return 0.0121 + 5.92e-05 * p - 2.61e-08 * p**2 + 5.11e-12 * p**3
        elif year=="2018":  return 0.0136 + 5.47e-05 * p - 2.3e-08  * p**2 + 4.66e-12 * p**3
    else:
        return 0

def get_highpt_mu_reso_smearing_sfs(year,p,eta):
    sigma_p = muon_momentum_reso_sigma(year, p, eta)
    nocorr = 1.0
    fiveperc = 1 + random.gauss(0, sigma_p * 0.3202)
    tenperc = 1 + random.gauss(0, sigma_p * 0.46)

    if abs(eta) < 1.2:
        if year=="2016APV": return( nocorr, tenperc ) #no additional smearing required for sf, 10% smearing required for syst calc
        elif year=="2016":  return( nocorr, tenperc ) #no additional smearing required for sf, 10% smearing required for syst calc
        elif year=="2017":  return( nocorr, tenperc ) #no additional smearing required for sf, 10% smearing required for syst calc
        elif year=="2018":  return( nocorr, tenperc ) #no additional smearing required for sf, 10% smearing required for syst calc
    elif 1.2 <= abs(eta) < 2.4:
        if year=="2016APV": return( nocorr,  tenperc ) #no additional smearing required for sf, 10% smearing required for syst calc
        elif year=="2016":  return( nocorr,  tenperc ) #no additional smearing required for sf, 10% smearing required for syst calc
        elif year=="2017":  return( fiveperc, 1 + random.gauss(0, fiveperc * 0.46) ) #5% smearing required for sf, 10% smearing required for syst calc
        elif year=="2018":  return( tenperc,  1 + random.gauss(0, tenperc * 0.46) ) #10% smearing required for sf, 10% smearing required for syst calc
    else:
        return( nocorr, nocorr )
        
def get_mu_scale_factor(year, categ, abseta, pt, scale_factor_type): #is actually p, not pt, for muon_reco
    if year=="2016APV": year="2016_preVFP"

    if categ=="reco":
        correction_file = "MuonSFs/ScaleFactors_Muon_highPt_RECO_"+year+"_schemaV2.json"
        correction_set  = CorrectionSet.from_file(correction_file)
        correction      = correction_set['NUM_GlobalMuons_DEN_TrackerMuonProbes']
    elif categ=="id":
        correction_file = "MuonSFs/ScaleFactors_Muon_highPt_IDISO_"+year+"_schemaV2.json"
        correction_set  = CorrectionSet.from_file(correction_file)
        correction      = correction_set['NUM_HighPtID_DEN_GlobalMuonProbes']
    elif categ=="iso":
        correction_file = "MuonSFs/ScaleFactors_Muon_highPt_IDISO_"+year+"_schemaV2.json"
        correction_set  = CorrectionSet.from_file(correction_file)
        correction      = correction_set['NUM_probe_LooseRelTkIso_DEN_HighPtProbes']
    elif categ=="hlt":
        correction_file = "MuonSFs/ScaleFactors_Muon_highPt_HLT_"+year+"_schemaV2.json"
        correction_set  = CorrectionSet.from_file(correction_file)
        correction      = correction_set['NUM_HLT_DEN_HighPtLooseRelIsoProbes']
    return correction.evaluate(abseta, pt, scale_factor_type)

def load_btagging_efficiency_histograms(year,mctype):
    if   int(mctype) in range(500,600): b_tagging_efficiency_file = "BTagEfficiencies/btag_efficiencies_2018_eta.root" #do a combo of all mass points (rn uses ?? (500 or 1000, probably 500))
    elif int(mctype) in range(600,700): b_tagging_efficiency_file = "BTagEfficiencies/btag_efficiencies_2018_etaprime.root" #do a combo of all mass points (rn uses M1000)
    elif int(mctype)==20:    b_tagging_efficiency_file = "BTagEfficiencies/btag_efficiencies_"+year+"_ttjets.root"
    elif int(mctype)==21:  b_tagging_efficiency_file = "BTagEfficiencies/btag_efficiencies_"+year+"_wjetstolnu.root"
    elif int(mctype)==22:  b_tagging_efficiency_file = "BTagEfficiencies/btag_efficiencies_"+year+"_dyjetstoll.root"
    
    myfile = ROOT.TFile.Open(b_tagging_efficiency_file, "READ")
    if not myfile or myfile.IsZombie():
        raise IOError(f"Could not open ROOT file: {root_file}")
    
    efficiency_b = myfile.Get("efficiency_b")
    efficiency_c = myfile.Get("efficiency_c")
    efficiency_light = myfile.Get("efficiency_light")

    if not efficiency_b or not efficiency_c or not efficiency_light:
        raise ValueError("Efficiency histograms not found in the ROOT file.")

    efficiency_b = efficiency_b.Clone()
    efficiency_c = efficiency_c.Clone()
    efficiency_light = efficiency_light.Clone()

    efficiency_b.SetDirectory(0)
    efficiency_c.SetDirectory(0)
    efficiency_light.SetDirectory(0)

    myfile.Close()

    return efficiency_b, efficiency_c, efficiency_light

def get_btagging_efficiency(pt, eta, flavour, efficiency_b, efficiency_c, efficiency_light):
    if flavour == 5:  # b-jets
        hist = efficiency_b
    elif flavour == 4:  # c-jets
        hist = efficiency_c
    elif flavour == 0:  # light-jets
        hist = efficiency_light
    else:
        raise ValueError("Invalid jet flavour. Use 5 for b, 4 for c, and 0 for light.")

    # Get the bin indices for pt and eta
    x_bin = hist.GetXaxis().FindBin(eta)
    y_bin = hist.GetYaxis().FindBin(pt)

    # Retrieve efficiency
    if 1 <= x_bin <= hist.GetNbinsX() and 1 <= y_bin <= hist.GetNbinsY():
        efficiency = hist.GetBinContent(x_bin, y_bin)
        return efficiency
    else:
        print(f"Warning: pt={pt}, eta={eta} out of histogram range!")
        return 0.0  # Return 0 efficiency if out of range
    
def check_goodtwoprong(twoprongs): #checks if there is goodtwoprong, and returns leading symiso 2p if it exists
    good_2ps=[twoprong for twoprong in twoprongs if(twoprong.pt > 20 and abs(twoprong.eta) < 1.3)]
    if not good_2ps:
        return False, None

    good_symiso_2ps = [twoprong for twoprong in twoprongs if(twoprong.isTight)]
    if not good_symiso_2ps:
        return True, None

    leading_symiso_2p = max(good_symiso_2ps, key=lambda twoprong: twoprong.pt)
    return True, leading_symiso_2p

def check_goodmuon(muons): #checks if there is goodmuon (w/o isolation cut), and returns leading iso good muon pt if it exists
    good_muons = [muon for muon in muons if ( muon.tunepRelPt*muon.pt > 52 and abs(muon.eta) < 2.4 and muon.highPtId == 2 )]
    if not good_muons:
        return False, None
    
    good_iso_muons = [muon for muon in good_muons if muon.tkRelIso < 0.1]
    if not good_iso_muons:
        return True, None
    
    leading_iso_muon = max(good_iso_muons, key=lambda muon: muon.tunepRelPt*muon.pt)
    return True, leading_iso_muon

def deltaR(obj1, obj2):
    """Calculate deltaR between two objects with .eta and .phi."""
    delta_eta = obj1.eta - obj2.eta
    delta_phi = math.atan2(math.sin(obj1.phi - obj2.phi), math.cos(obj1.phi - obj2.phi))  # Correct for phi wrapping
    return math.sqrt(delta_eta**2 + delta_phi**2)

def check_goodjets(jets,leading_iso_muon,leading_symiso_2p): #check if there are >=3 goodjets, and returns collection of excellentjets if there are >=3
    if len(jets) < 3:
        return False, None

    goodjets = [jet for jet in jets if jet.pt > 30 and abs(jet.eta) < 2.5]
    if len(goodjets) < 3:
        return False, None

    excellentjets = [jet for jet in goodjets
                     if ( ((leading_iso_muon is None) or deltaR(jet, leading_iso_muon) >= 0.3)
                          and ((leading_symiso_2p is None) or deltaR(jet, leading_symiso_2p) >= 0.3) )]
    if len(excellentjets) < 3:
        return True, None

    return True, excellentjets 
    #19-11-2024: clean goodjets of the leading iso muon if it exists, and the leading isosym 2p if it exists

class attoAOD_ttw_mu(Module):
    def __init__(self, year="2018", mctype="0", attoVersion=version, run=None): 
        self.year = year
        self.mctype = mctype
        self.run = run
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
        if int(self.mctype) in range(500,700) or int(self.mctype) in range(20,30):
            self.out.branch("SF_MuonReco_nominal","F")
            self.out.branch("SF_MuonReco_up", "F")
            self.out.branch("SF_MuonReco_down","F")
            self.out.branch("SF_MuonId_nominal","F")
            self.out.branch("SF_MuonId_up", "F")
            self.out.branch("SF_MuonId_down","F")
            self.out.branch("SF_MuonIso_nominal","F")
            self.out.branch("SF_MuonIso_up", "F")
            self.out.branch("SF_MuonIso_down","F")
            self.out.branch("SF_MuonHlt_nominal","F")
            self.out.branch("SF_MuonHlt_up", "F")
            self.out.branch("SF_MuonHlt_down","F")
            self.out.branch("SF_BTagging_nominal","F")
            self.out.branch("SF_BTagging_bc_correlated_up","F")
            self.out.branch("SF_BTagging_bc_correlated_down","F")
            self.out.branch("SF_BTagging_light_correlated_up","F")
            self.out.branch("SF_BTagging_light_correlated_down","F")
            self.out.branch("SF_BTagging_bc_uncorrelated_up","F")
            self.out.branch("SF_BTagging_bc_uncorrelated_down","F")
            self.out.branch("SF_BTagging_light_uncorrelated_up","F")
            self.out.branch("SF_BTagging_light_uncorrelated_down","F")
            self.out.branch("SF_Pileup_nominal","F")
            self.out.branch("SF_Pileup_up", "F")
            self.out.branch("SF_Pileup_down","F")
            self.out.branch("Corr_MuonReso_nominal","F")
            self.out.branch("Corr_MuonReso_syst", "F")
        self.out.branch("Corr_MuonRoc_nominal","F")
        self.out.branch("Corr_MuonRoc_up", "F")
        self.out.branch("Corr_MuonRoc_down","F")
        self.out.branch("year","I")
        self.out.branch("mcType","I")
        self.out.branch("passTrigger","O")
        self.out.branch("attoVersion","I")

        #for jet corrs
        self.out.branch("nJetCorr",             "I")
        self.out.branch("JetCorr_pt",           "F", lenVar="nJetCorr")
        self.out.branch("JetCorr_corrval",      "F", lenVar="nJetCorr")
        self.out.branch("JetCorr_uncorrpt",     "F", lenVar="nJetCorr")
        if int(self.mctype) !=0:
            self.out.branch("Corr_JetCorr_Absolute_correlated",             "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_FlavorQCD_correlated",            "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_BBEC1_correlated",                "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_EC2_correlated",                  "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_HF_correlated",                   "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_RelativeBal_correlated",          "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_Absolute_uncorrelated",           "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_BBEC1_uncorrelated",              "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_EC2_uncorrelated",                "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_HF_uncorrelated",                 "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_RelativeSample_uncorrelated",     "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_Total",                           "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_JER",                                "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_JER_uncorrelated_up",                "F", lenVar="nJetCorr")
            self.out.branch("Corr_JetCorr_JER_uncorrelated_down",              "F", lenVar="nJetCorr")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        flags = Object(event, "Flag")
        trigger = Object(event, "HLT")
        jets = Collection(event, "Jet")
        twoprongs = Collection(event, "TwoProng")
        muons = Collection(event, "Muon")
        if self.mctype!="0":
            genparts = Collection(event, "GenPart")
            genjets = Collection(event, "GenJet")
        
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
            if (self.mygetattr(trigger, 'Mu50', False)) or (self.mygetattr(trigger, 'TkMu100', False)) or (self.mygetattr(trigger, 'OldMu100', False)):
                self.out.fillBranch("passTrigger", True)
            else: 
                self.out.fillBranch("passTrigger", False)
                if selections: return False
        elif (self.year=="2016" or self.year=="2016APV"):
            if (self.mygetattr(trigger, 'Mu50', False)) or (self.mygetattr(trigger, 'TkMu50', False)):
                self.out.fillBranch("passTrigger", True)
            else: 
                self.out.fillBranch("passTrigger", False)
                if selections: return False

        #event selection
        good_muon_exists, leading_iso_muon = check_goodmuon(muons)
        good_2p_exists, leading_symiso_2p = check_goodtwoprong(twoprongs)
        enough_goodjets_exist, cleaned_jets = check_goodjets(jets,leading_iso_muon,leading_symiso_2p)

        if selections:
            if not ( good_2p_exists and good_muon_exists and enough_goodjets_exist ): return False

        jetcorr_pt=[]
        jetcorr_corrval=[]
        jetcorr_uncorrpt=[]
        jetcorr_Absolute_correlated             =[]
        jetcorr_FlavorQCD_correlated            =[]
        jetcorr_BBEC1_correlated                =[]
        jetcorr_EC2_correlated                  =[]
        jetcorr_HF_correlated                   =[]
        jetcorr_RelativeBal_correlated          =[]
        jetcorr_Absolute_uncorrelated           =[]
        jetcorr_BBEC1_uncorrelated              =[]
        jetcorr_EC2_uncorrelated                =[]
        jetcorr_HF_uncorrelated                 =[]
        jetcorr_RelativeSample_uncorrelated     =[]
        jetcorr_Total                           =[]
        jetcorr_JER                             =[]
        jetcorr_JER_uncorrelated_up             =[]
        jetcorr_JER_uncorrelated_down           =[]
        
        for jet in jets:
            genjetPt=-999
            if self.mctype!="0" and jet.genJetIdx>=0: #matched mc jet
                try:
                    genjetPt = genjets[(jet.genJetIdx)].pt
                except IndexError:
                    genjetPt=-999

            if (not enough_goodjets_exist) or (jet.pt < 30) or (abs(jet.eta) > 2.5):
                jetcorr_pt.append(jet.pt)
                jetcorr_corrval.append(1)
                jetcorr_uncorrpt.append(jet.pt)
                if self.mctype!="0":
                    jetcorr_Absolute_correlated.append(0)
                    jetcorr_FlavorQCD_correlated.append(0)
                    jetcorr_BBEC1_correlated.append(0)
                    jetcorr_EC2_correlated.append(0)
                    jetcorr_HF_correlated.append(0)
                    jetcorr_RelativeBal_correlated.append(0)
                    jetcorr_Absolute_uncorrelated.append(0)
                    jetcorr_BBEC1_uncorrelated.append(0)
                    jetcorr_EC2_uncorrelated.append(0)
                    jetcorr_HF_uncorrelated.append(0)
                    jetcorr_RelativeSample_uncorrelated.append(0)
                    jetcorr_Total.append(0)
                    jetcorr_JER.append(1)
                    jetcorr_JER_uncorrelated_up.append(1)
                    jetcorr_JER_uncorrelated_down.append(1)
            else:
                if self.mctype=="0":
                    corrval,pt_corr=correct_jet(self.year, self.mctype, self.run, jet.rawFactor, jet.pt, jet.eta, jet.mass, jet.area, event.fixedGridRhoFastjetAll, genjetPt)
                    jetcorr_pt.append(pt_corr)
                    jetcorr_corrval.append(corrval)
                    jetcorr_uncorrpt.append(jet.pt)
                else:
                    corrval,pt_corr, reso_sf,reso_sf_up,reso_sf_down, correl_uncert_vals,uncorrel_uncert_vals = correct_jet(self.year, self.mctype, self.run, jet.rawFactor, jet.pt, jet.eta, jet.mass, jet.area, event.fixedGridRhoFastjetAll, genjetPt)
                    jetcorr_pt.append(pt_corr)
                    jetcorr_corrval.append(corrval)
                    jetcorr_uncorrpt.append(jet.pt)
                    jetcorr_Absolute_correlated.append(correl_uncert_vals[0])
                    jetcorr_FlavorQCD_correlated.append(correl_uncert_vals[1])
                    jetcorr_BBEC1_correlated.append(correl_uncert_vals[2])
                    jetcorr_EC2_correlated.append(correl_uncert_vals[3])
                    jetcorr_HF_correlated.append(correl_uncert_vals[4])
                    jetcorr_RelativeBal_correlated.append(correl_uncert_vals[5])
                    jetcorr_Total.append(correl_uncert_vals[6])
                    jetcorr_Absolute_uncorrelated.append(uncorrel_uncert_vals[0])
                    jetcorr_BBEC1_uncorrelated.append(uncorrel_uncert_vals[1])
                    jetcorr_EC2_uncorrelated.append(uncorrel_uncert_vals[2])
                    jetcorr_HF_uncorrelated.append(uncorrel_uncert_vals[3])
                    jetcorr_RelativeSample_uncorrelated.append(uncorrel_uncert_vals[4])
                    jetcorr_JER.append(reso_sf)
                    jetcorr_JER_uncorrelated_up.append(reso_sf_up)
                    jetcorr_JER_uncorrelated_down.append(reso_sf_down)
                
        #Construct b-tagging related event weights
        if ( int(self.mctype) in range(500,700) or int(self.mctype) in range(20,30))  and cleaned_jets is not None:
            efficiency_b, efficiency_c, efficiency_light = load_btagging_efficiency_histograms(self.year,self.mctype)
            btagged=[]
            efficiency=[]
            hadronflav=[]
            if self.year=='2018':
                btagMwp=0.2783
                correction_set_btv = CorrectionSet.from_file('POG/BTV/2018_UL/btagging.json.gz')
            elif self.year=='2017':
                btagMwp=0.3040
                correction_set_btv = CorrectionSet.from_file('POG/BTV/2017_UL/btagging.json.gz')
            elif self.year=='2016':
                btagMwp=0.2489
                correction_set_btv = CorrectionSet.from_file('POG/BTV/2016_UL/btagging.json.gz')
            elif self.year=='2016APV':
                btagMwp=0.2598
                correction_set_btv = CorrectionSet.from_file('POG/BTV/2016APV_UL/btagging.json.gz')
            central=[]
            up=[]
            down=[]
            upcorrelated=[]
            downcorrelated=[]
            for jet in cleaned_jets:
                btagged.append(jet.btagDeepFlavB > btagMwp)
                efficiency.append( get_btagging_efficiency(jet.pt, jet.eta, jet.hadronFlavour, efficiency_b, efficiency_c, efficiency_light) )
                hadronflav.append( jet.hadronFlavour )
                if jet.hadronFlavour==0:
                    central.append(       correction_set_btv["deepJet_incl"].evaluate("central",        "M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                    up.append(            correction_set_btv["deepJet_incl"].evaluate("up",             "M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                    down.append(          correction_set_btv["deepJet_incl"].evaluate("down",           "M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                    upcorrelated.append(  correction_set_btv["deepJet_incl"].evaluate("up_correlated",  "M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                    downcorrelated.append(correction_set_btv["deepJet_incl"].evaluate("down_correlated","M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                else:
                    central.append(       correction_set_btv["deepJet_comb"].evaluate("central",        "M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                    up.append(            correction_set_btv["deepJet_comb"].evaluate("up",             "M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                    down.append(          correction_set_btv["deepJet_comb"].evaluate("down",           "M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                    upcorrelated.append(  correction_set_btv["deepJet_comb"].evaluate("up_correlated",  "M",jet.hadronFlavour,abs(jet.eta),jet.pt))
                    downcorrelated.append(correction_set_btv["deepJet_comb"].evaluate("down_correlated","M",jet.hadronFlavour,abs(jet.eta),jet.pt))

            nominal=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                           (jet_idx for jet_idx, is_btagged in enumerate(btagged) if is_btagged),1.0) *\
                reduce(lambda prod, jet_idx:      prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if not is_btagged),1.0)

            bc_correlated_up=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                                (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]==0))    ,1.0) *\
                reduce(lambda prod, jet_idx:              prod * upcorrelated[jet_idx],                                                           (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]!=0))    ,1.0) *\
                reduce(lambda prod, jet_idx:              prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])),      (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]==0)),1.0) *\
                reduce(lambda prod, jet_idx:              prod * ((1 - upcorrelated[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]!=0)),1.0)

            bc_correlated_down=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                                  (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]==0))    ,1.0) *\
                reduce(lambda prod, jet_idx:                prod * downcorrelated[jet_idx],                                                           (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]!=0))    ,1.0) *\
                reduce(lambda prod, jet_idx:                prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])),        (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]==0)),1.0) *\
                reduce(lambda prod, jet_idx:                prod * ((1 - downcorrelated[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]!=0)),1.0)

            light_correlated_up=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                             (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]!=0))    ,1.0) *\
                reduce(lambda prod, jet_idx:              prod * upcorrelated[jet_idx],                                                           (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]==0))    ,1.0) *\
                reduce(lambda prod, jet_idx:              prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])),      (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]!=0)),1.0) *\
                reduce(lambda prod, jet_idx:              prod * ((1 - upcorrelated[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]==0)),1.0)

            light_correlated_down=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                               (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]!=0))    ,1.0) *\
                reduce(lambda prod, jet_idx:                prod * downcorrelated[jet_idx],                                                           (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]==0))    ,1.0) *\
                reduce(lambda prod, jet_idx:                prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])),        (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]!=0)),1.0) *\
                reduce(lambda prod, jet_idx:                prod * ((1 - downcorrelated[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]==0)),1.0)

            bc_uncorrelated_up=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                         (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]==0))    ,1.0) *\
                reduce(lambda prod, jet_idx:              prod * up[jet_idx],                                                                (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]!=0))    ,1.0) *\
                reduce(lambda prod, jet_idx:              prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]==0)),1.0) *\
                reduce(lambda prod, jet_idx:              prod * ((1 - up[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])),      (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]!=0)),1.0)

            bc_uncorrelated_down=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                         (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]==0))    ,1.0) *\
                reduce(lambda prod, jet_idx:                prod * down[jet_idx],                                                              (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]!=0))    ,1.0) *\
                reduce(lambda prod, jet_idx:                prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]==0)),1.0) *\
                reduce(lambda prod, jet_idx:                prod * ((1 - down[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])),    (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]!=0)),1.0)

            light_uncorrelated_up=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                      (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]!=0))    ,1.0) *\
                reduce(lambda prod, jet_idx:              prod * up[jet_idx],                                                                (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]==0))    ,1.0) *\
                reduce(lambda prod, jet_idx:              prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]!=0)),1.0) *\
                reduce(lambda prod, jet_idx:              prod * ((1 - up[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])),      (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]==0)),1.0)

            light_uncorrelated_down=reduce(lambda prod, jet_idx: prod * central[jet_idx],                                                      (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]!=0))    ,1.0) *\
                reduce(lambda prod, jet_idx:                prod * down[jet_idx],                                                              (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (is_btagged and hadronflav[jet_idx]==0))    ,1.0) *\
                reduce(lambda prod, jet_idx:                prod * ((1 - central[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])), (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]!=0)),1.0) *\
                reduce(lambda prod, jet_idx:                prod * ((1 - down[jet_idx] * efficiency[jet_idx]) / (1 - efficiency[jet_idx])),    (jet_idx for jet_idx, is_btagged in enumerate(btagged) if (not is_btagged and hadronflav[jet_idx]==0)),1.0)

            #print(nominal,bc_correlated_up,bc_correlated_down,light_correlated_up,light_correlated_down,bc_uncorrelated_up,bc_uncorrelated_down,light_uncorrelated_up,light_uncorrelated_down)
            
        #fill branches
        if self.year != "2016APV": self.out.fillBranch("year", int(self.year))
        else: self.out.fillBranch("year", int(20160))
        self.out.fillBranch("mcType", int(self.mctype))
        self.out.fillBranch("attoVersion", int(self.attoVersion))

        self.out.fillBranch("nJetCorr",    len(jetcorr_pt))
        self.out.fillBranch("JetCorr_pt",       jetcorr_pt)
        self.out.fillBranch("JetCorr_corrval",  jetcorr_corrval)
        self.out.fillBranch("JetCorr_uncorrpt", jetcorr_uncorrpt)
        if int(self.mctype) !=0:
            self.out.fillBranch("Corr_JetCorr_Absolute_correlated",        jetcorr_Absolute_correlated)
            self.out.fillBranch("Corr_JetCorr_FlavorQCD_correlated",       jetcorr_FlavorQCD_correlated)
            self.out.fillBranch("Corr_JetCorr_BBEC1_correlated",           jetcorr_BBEC1_correlated)
            self.out.fillBranch("Corr_JetCorr_EC2_correlated",             jetcorr_EC2_correlated)
            self.out.fillBranch("Corr_JetCorr_HF_correlated",              jetcorr_HF_correlated)
            self.out.fillBranch("Corr_JetCorr_RelativeBal_correlated",     jetcorr_RelativeBal_correlated)
            self.out.fillBranch("Corr_JetCorr_Total",                      jetcorr_Total)
            self.out.fillBranch("Corr_JetCorr_Absolute_uncorrelated",      jetcorr_Absolute_uncorrelated)
            self.out.fillBranch("Corr_JetCorr_BBEC1_uncorrelated",         jetcorr_BBEC1_uncorrelated)
            self.out.fillBranch("Corr_JetCorr_EC2_uncorrelated",           jetcorr_EC2_uncorrelated)
            self.out.fillBranch("Corr_JetCorr_HF_uncorrelated",            jetcorr_HF_uncorrelated)
            self.out.fillBranch("Corr_JetCorr_RelativeSample_uncorrelated",jetcorr_RelativeSample_uncorrelated)
            self.out.fillBranch("Corr_JetCorr_JER",                        jetcorr_JER)
            self.out.fillBranch("Corr_JetCorr_JER_uncorrelated_up",        jetcorr_JER_uncorrelated_up)
            self.out.fillBranch("Corr_JetCorr_JER_uncorrelated_down",      jetcorr_JER_uncorrelated_down)

        if leading_iso_muon:
            genPt=-999
            if self.mctype!="0" and leading_iso_muon.genPartFlav !=0: #matched mc muon
                genPt = genparts[leading_iso_muon.genPartIdx].pt
            roc_sf, roc_err= get_rochester_mu_corrections(self.year, self.mctype, leading_iso_muon, genPt)
            self.out.fillBranch("Corr_MuonRoc_nominal", roc_sf)
            self.out.fillBranch("Corr_MuonRoc_up",      roc_sf + roc_err)
            self.out.fillBranch("Corr_MuonRoc_down",    roc_sf - roc_err)
        else:
            self.out.fillBranch("Corr_MuonRoc_nominal", 1.0)
            self.out.fillBranch("Corr_MuonRoc_up",      1.0)
            self.out.fillBranch("Corr_MuonRoc_down",    1.0)
                
        if int(self.mctype) in range(500,700) or int(self.mctype) in range(20,30): #only for MC
            pu_weight, pu_weight_up, pu_weight_down = 1.0, 1.0, 1.0
            pu_weight, pu_weight_up, pu_weight_down = get_pileup_weight(self.year, event.Pileup_nTrueInt)
            self.out.fillBranch("SF_Pileup_nominal", pu_weight)
            self.out.fillBranch("SF_Pileup_up",      pu_weight_up)
            self.out.fillBranch("SF_Pileup_down",    pu_weight_down)
        
            if leading_iso_muon:
                leading_iso_muon_vec = ROOT.TLorentzVector()
                leading_iso_muon_vec.SetPtEtaPhiM(leading_iso_muon.tunepRelPt*leading_iso_muon.pt, leading_iso_muon.eta, leading_iso_muon.phi, leading_iso_muon.mass)
                leading_iso_muon_p=leading_iso_muon_vec.P()
                self.out.fillBranch("SF_MuonReco_nominal", get_mu_scale_factor(self.year, "reco", abs(leading_iso_muon.eta), leading_iso_muon_p, "nominal"))
                self.out.fillBranch("SF_MuonReco_up",      get_mu_scale_factor(self.year, "reco", abs(leading_iso_muon.eta), leading_iso_muon_p, "systup"))
                self.out.fillBranch("SF_MuonReco_down",    get_mu_scale_factor(self.year, "reco", abs(leading_iso_muon.eta), leading_iso_muon_p, "systdown"))
                self.out.fillBranch("SF_MuonId_nominal",   get_mu_scale_factor(self.year, "id", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "nominal"))
                self.out.fillBranch("SF_MuonId_up",        get_mu_scale_factor(self.year, "id", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "systup"))
                self.out.fillBranch("SF_MuonId_down",      get_mu_scale_factor(self.year, "id", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "systdown"))
                self.out.fillBranch("SF_MuonIso_nominal",  get_mu_scale_factor(self.year, "iso", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "nominal"))
                self.out.fillBranch("SF_MuonIso_up",       get_mu_scale_factor(self.year, "iso", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "systup"))
                self.out.fillBranch("SF_MuonIso_down",     get_mu_scale_factor(self.year, "iso", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "systdown"))
                self.out.fillBranch("SF_MuonHlt_nominal",  get_mu_scale_factor(self.year, "hlt", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "nominal"))
                self.out.fillBranch("SF_MuonHlt_up",       get_mu_scale_factor(self.year, "hlt", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "systup"))
                self.out.fillBranch("SF_MuonHlt_down",     get_mu_scale_factor(self.year, "hlt", abs(leading_iso_muon.eta), leading_iso_muon.tunepRelPt*leading_iso_muon.pt, "systdown"))
                muon_reso_sf, muon_reso_syst_sf = get_highpt_mu_reso_smearing_sfs(self.year, leading_iso_muon_p, leading_iso_muon.eta)
                self.out.fillBranch("Corr_MuonReso_nominal", muon_reso_sf)
                self.out.fillBranch("Corr_MuonReso_syst",    muon_reso_syst_sf)
            else:
                self.out.fillBranch("SF_MuonReco_nominal", -99.0)
                self.out.fillBranch("SF_MuonReco_up",      -99.0)
                self.out.fillBranch("SF_MuonReco_down",    -99.0)
                self.out.fillBranch("SF_MuonId_nominal",   -99.0)
                self.out.fillBranch("SF_MuonId_up",        -99.0)
                self.out.fillBranch("SF_MuonId_down",      -99.0)
                self.out.fillBranch("SF_MuonIso_nominal",  -99.0)
                self.out.fillBranch("SF_MuonIso_up",       -99.0)
                self.out.fillBranch("SF_MuonIso_down",     -99.0)
                self.out.fillBranch("SF_MuonHlt_nominal",  -99.0)
                self.out.fillBranch("SF_MuonHlt_up",       -99.0)
                self.out.fillBranch("SF_MuonHlt_down",     -99.0)
                self.out.fillBranch("Corr_MuonReso_nominal", 1.0)
                self.out.fillBranch("Corr_MuonReso_syst",    1.0)
            if cleaned_jets:
                self.out.fillBranch("SF_BTagging_nominal",                 nominal)
                self.out.fillBranch("SF_BTagging_bc_correlated_up",        bc_correlated_up)
                self.out.fillBranch("SF_BTagging_bc_correlated_down",      bc_correlated_down)
                self.out.fillBranch("SF_BTagging_light_correlated_up",     light_correlated_up)
                self.out.fillBranch("SF_BTagging_light_correlated_down",   light_correlated_down)
                self.out.fillBranch("SF_BTagging_bc_uncorrelated_up",      bc_uncorrelated_up)
                self.out.fillBranch("SF_BTagging_bc_uncorrelated_down",    bc_uncorrelated_down)
                self.out.fillBranch("SF_BTagging_light_uncorrelated_up",   light_uncorrelated_up)
                self.out.fillBranch("SF_BTagging_light_uncorrelated_down", light_uncorrelated_down)
            else:
                self.out.fillBranch("SF_BTagging_nominal",                 -99.0)
                self.out.fillBranch("SF_BTagging_bc_correlated_up",        -99.0)
                self.out.fillBranch("SF_BTagging_bc_correlated_down",      -99.0)
                self.out.fillBranch("SF_BTagging_light_correlated_up",     -99.0)
                self.out.fillBranch("SF_BTagging_light_correlated_down",   -99.0)
                self.out.fillBranch("SF_BTagging_bc_uncorrelated_up",      -99.0)
                self.out.fillBranch("SF_BTagging_bc_uncorrelated_down",    -99.0)
                self.out.fillBranch("SF_BTagging_light_uncorrelated_up",   -99.0)
                self.out.fillBranch("SF_BTagging_light_uncorrelated_down", -99.0)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
mu_2018_singlemuonA =      lambda: attoAOD_ttw_mu(run="A")
mu_2018_singlemuonB =      lambda: attoAOD_ttw_mu(run="B")
mu_2018_singlemuonC =      lambda: attoAOD_ttw_mu(run="C")
mu_2018_singlemuonD =      lambda: attoAOD_ttw_mu(run="D")
mu_2018_ttjets =           lambda: attoAOD_ttw_mu(mctype="20")
mu_2018_wjetstolnu =       lambda: attoAOD_ttw_mu(mctype="21")
mu_2018_dyjetstoll =       lambda: attoAOD_ttw_mu(mctype="22")
mu_2018_eta_M500  =        lambda: attoAOD_ttw_mu(mctype="500")
mu_2018_eta_M750  =        lambda: attoAOD_ttw_mu(mctype="501")
mu_2018_eta_M850  =        lambda: attoAOD_ttw_mu(mctype="502")
mu_2018_eta_M1000 =        lambda: attoAOD_ttw_mu(mctype="503")
mu_2018_eta_M1500 =        lambda: attoAOD_ttw_mu(mctype="504")
mu_2018_eta_M2000 =        lambda: attoAOD_ttw_mu(mctype="505")
mu_2018_eta_M2500 =        lambda: attoAOD_ttw_mu(mctype="506")
mu_2018_eta_M3000 =        lambda: attoAOD_ttw_mu(mctype="507")
mu_2018_eta_M4000 =        lambda: attoAOD_ttw_mu(mctype="508")
mu_2018_etaprime_M850  =   lambda: attoAOD_ttw_mu(mctype="602")
mu_2018_etaprime_M1000 =   lambda: attoAOD_ttw_mu(mctype="603")
mu_2018_etaprime_M1500 =   lambda: attoAOD_ttw_mu(mctype="604")
mu_2018_etaprime_M2000 =   lambda: attoAOD_ttw_mu(mctype="605")
mu_2018_etaprime_M2500 =   lambda: attoAOD_ttw_mu(mctype="606")
mu_2018_etaprime_M3000 =   lambda: attoAOD_ttw_mu(mctype="607")
mu_2018_etaprime_M4000 =   lambda: attoAOD_ttw_mu(mctype="608")

mu_2017_singlemuonB =      lambda: attoAOD_ttw_mu(year="2017",run="B")
mu_2017_singlemuonC =      lambda: attoAOD_ttw_mu(year="2017",run="C")
mu_2017_singlemuonD =      lambda: attoAOD_ttw_mu(year="2017",run="D")
mu_2017_singlemuonE =      lambda: attoAOD_ttw_mu(year="2017",run="E")
mu_2017_singlemuonF =      lambda: attoAOD_ttw_mu(year="2017",run="F")
mu_2017_ttjets =           lambda: attoAOD_ttw_mu(year="2017",mctype="20")
mu_2017_wjetstolnu =       lambda: attoAOD_ttw_mu(year="2017",mctype="21")
mu_2017_dyjetstoll =       lambda: attoAOD_ttw_mu(year="2017",mctype="22")
mu_2017_eta_M500  =        lambda: attoAOD_ttw_mu(year="2017",mctype="500")
mu_2017_eta_M750  =        lambda: attoAOD_ttw_mu(year="2017",mctype="501")
mu_2017_eta_M850  =        lambda: attoAOD_ttw_mu(year="2017",mctype="502")
mu_2017_eta_M1000 =        lambda: attoAOD_ttw_mu(year="2017",mctype="503")
mu_2017_eta_M1500 =        lambda: attoAOD_ttw_mu(year="2017",mctype="504")
mu_2017_eta_M2000 =        lambda: attoAOD_ttw_mu(year="2017",mctype="505")
mu_2017_eta_M2500 =        lambda: attoAOD_ttw_mu(year="2017",mctype="506")
mu_2017_eta_M3000 =        lambda: attoAOD_ttw_mu(year="2017",mctype="507")
mu_2017_eta_M4000 =        lambda: attoAOD_ttw_mu(year="2017",mctype="508")
mu_2017_etaprime_M850  =   lambda: attoAOD_ttw_mu(year="2017",mctype="602")
mu_2017_etaprime_M1000 =   lambda: attoAOD_ttw_mu(year="2017",mctype="603")
mu_2017_etaprime_M1500 =   lambda: attoAOD_ttw_mu(year="2017",mctype="604")
mu_2017_etaprime_M2000 =   lambda: attoAOD_ttw_mu(year="2017",mctype="605")
mu_2017_etaprime_M2500 =   lambda: attoAOD_ttw_mu(year="2017",mctype="606")
mu_2017_etaprime_M3000 =   lambda: attoAOD_ttw_mu(year="2017",mctype="607")
mu_2017_etaprime_M4000 =   lambda: attoAOD_ttw_mu(year="2017",mctype="608")

mu_2016_singlemuonF2 =     lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonG =      lambda: attoAOD_ttw_mu(year="2016")
mu_2016_singlemuonH =      lambda: attoAOD_ttw_mu(year="2016")
mu_2016_ttjets =           lambda: attoAOD_ttw_mu(year="2016",mctype="20")
mu_2016_wjetstolnu =       lambda: attoAOD_ttw_mu(year="2016",mctype="21")
mu_2016_dyjetstoll =       lambda: attoAOD_ttw_mu(year="2016",mctype="22")
mu_2016_eta_M500  =        lambda: attoAOD_ttw_mu(year="2016",mctype="500")
mu_2016_eta_M750  =        lambda: attoAOD_ttw_mu(year="2016",mctype="501")
mu_2016_eta_M850  =        lambda: attoAOD_ttw_mu(year="2016",mctype="502")
mu_2016_eta_M1000 =        lambda: attoAOD_ttw_mu(year="2016",mctype="503")
mu_2016_eta_M1500 =        lambda: attoAOD_ttw_mu(year="2016",mctype="504")
mu_2016_eta_M2000 =        lambda: attoAOD_ttw_mu(year="2016",mctype="505")
mu_2016_eta_M2500 =        lambda: attoAOD_ttw_mu(year="2016",mctype="506")
mu_2016_eta_M3000 =        lambda: attoAOD_ttw_mu(year="2016",mctype="507")
mu_2016_eta_M4000 =        lambda: attoAOD_ttw_mu(year="2016",mctype="508")
mu_2016_etaprime_M850  =   lambda: attoAOD_ttw_mu(year="2016",mctype="602")
mu_2016_etaprime_M1000 =   lambda: attoAOD_ttw_mu(year="2016",mctype="603")
mu_2016_etaprime_M1500 =   lambda: attoAOD_ttw_mu(year="2016",mctype="604")
mu_2016_etaprime_M2000 =   lambda: attoAOD_ttw_mu(year="2016",mctype="605")
mu_2016_etaprime_M2500 =   lambda: attoAOD_ttw_mu(year="2016",mctype="606")
mu_2016_etaprime_M3000 =   lambda: attoAOD_ttw_mu(year="2016",mctype="607")
mu_2016_etaprime_M4000 =   lambda: attoAOD_ttw_mu(year="2016",mctype="608")

mu_2016APV_singlemuonB2 =     lambda: attoAOD_ttw_mu(year="2016APV",run="BCD")
mu_2016APV_singlemuonC =      lambda: attoAOD_ttw_mu(year="2016APV",run="BCD")
mu_2016APV_singlemuonD =      lambda: attoAOD_ttw_mu(year="2016APV",run="BCD")
mu_2016APV_singlemuonE =      lambda: attoAOD_ttw_mu(year="2016APV",run="EF")
mu_2016APV_singlemuonF =      lambda: attoAOD_ttw_mu(year="2016APV",run="EF")
mu_2016APV_ttjets =           lambda: attoAOD_ttw_mu(year="2016APV",mctype="20")
mu_2016APV_wjetstolnu =       lambda: attoAOD_ttw_mu(year="2016APV",mctype="21")
mu_2016APV_dyjetstoll =       lambda: attoAOD_ttw_mu(year="2016APV",mctype="22")
mu_2016APV_eta_M500  =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="500")
mu_2016APV_eta_M750  =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="501")
mu_2016APV_eta_M850  =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="502")
mu_2016APV_eta_M1000 =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="503")
mu_2016APV_eta_M1500 =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="504")
mu_2016APV_eta_M2000 =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="505")
mu_2016APV_eta_M2500 =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="506")
mu_2016APV_eta_M3000 =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="507")
mu_2016APV_eta_M4000 =        lambda: attoAOD_ttw_mu(year="2016APV",mctype="508")
mu_2016APV_etaprime_M850  =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="602")
mu_2016APV_etaprime_M1000 =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="603")
mu_2016APV_etaprime_M1500 =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="604")
mu_2016APV_etaprime_M2000 =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="605")
mu_2016APV_etaprime_M2500 =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="606")
mu_2016APV_etaprime_M3000 =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="607")
mu_2016APV_etaprime_M4000 =   lambda: attoAOD_ttw_mu(year="2016APV",mctype="608")
