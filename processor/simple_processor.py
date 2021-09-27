#!/usr/bin/env python
import coffea
import numpy as np
import awkward as ak
from coffea import hist, processor
import topcoffea.modules.objects as top


#def isTightElec(pt, eta, dxy, dz, miniIso, sip3D, mvaTTH, elecMVA, lostHits, convVeto, tightCharge, sieie, hoe, eInvMinusPInv, minpt=15.0):
#  maskPOGMVA = ((pt<10)&(abs(eta)<0.8)&(elecMVA>-0.13))|((pt<10)&(abs(eta)>0.8)&(abs(eta)<1.44)&(elecMVA>-0.32))|((pt<10)&(abs(eta)>1.44)&(elecMVA>-0.08))|\
#               ((pt>10)&(abs(eta)<0.8)&(elecMVA>-0.86))|((pt>10)&(abs(eta)>0.8)&(abs(eta)<1.44)&(elecMVA>-0.81))|((pt>10)&(abs(eta)>1.44)&(elecMVA>-0.72))
#  maskSieie  = ((abs(eta)<1.479)&(sieie<0.011))|((abs(eta)>1.479)&(sieie<0.030))
#  maskhoe    = ((abs(eta)<1.479)&(hoe<0.10))|((abs(eta)>1.479)&(hoe<0.07))
#  mask = (pt>minpt)&(abs(eta)<2.5)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(sip3D<8)&(lostHits<=1)&\
#         (convVeto)&(maskSieie)&(maskPOGMVA)&(eInvMinusPInv>-0.04)&(maskhoe)&(miniIso<0.25)&(mvaTTH>0.90)&(tightCharge==2)
#  return mask
  
#def isTightMuon(pt, eta, dxy, dz, miniIso, sip3D, mvaTTH, mediumPrompt, tightCharge, looseId, minpt=10.0):
#  mask = (pt>minpt)&(abs(eta)<2.5)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(sip3D<8)&(looseId)&(miniIso<0.25)&(mvaTTH>0.90)&(tightCharge==2)&(mediumPrompt)
#  return mask

def isPresElec(pt, eta, dxy, dz, miniIso):
    pt    = (pt       > 32)
    eta   = (abs(eta) < 2.4)
    dxy   = (abs(dxy) < 0.05)
    dz    = (abs(dz)  < 0.1)
    iso   = (miniIso  < 0.4)
    return (pt & eta & dxy & dz & iso)

def isPresMuon(dxy, dz, eta, pt):
    pt_mask    = (pt         > 27)
    eta_mask   = (abs(eta)   < 2.4)
    dxy_mask   = (abs(dxy)   < 0.2)
    dz_mask    = (abs(dz)    < 0.5)
    return (pt_mask & eta_mask & dxy_mask & dz_mask)

def isLooseElec(miniPFRelIso_all,sip3d,lostHits):
    return (miniPFRelIso_all<0.4) & (sip3d<8) & (lostHits<=1)

def isLooseMuon(miniPFRelIso_all,sip3d):
    return (miniPFRelIso_all<0.4) & (sip3d<8)

def isClean(obj_A, obj_B, drmin=0.4):
    objB_near, objB_DR = obj_A.nearest(obj_B, return_metric=True)
    mask = ak.fill_none(objB_DR > drmin, True)
    return (mask)
	

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self):

        # Create the histograms
        # In general, histograms depend on 'sample', 'channel' (final state) and 'cut' (level of selection)
        self._accumulator = processor.dict_accumulator({
            "events": hist.Hist(
                  "Events",
                  hist.Cat("sample", "A"),
                  hist.Bin("pdgID", "B", 1200, -600, 600 )),
            "electron_pt": hist.Hist(
                  "Events",
                  hist.Cat("sample", "Electron p_t"),
                  hist.Bin("electron_pt", "Electron p_t", 60, 0, 600)),
            "muon_pt": hist.Hist(
                  "Events",
                  hist.Cat("sample", "Muon p_t"),
                  hist.Bin("muon_pt", "Muon p_t", 60, 0, 600)),
            "jet_pt": hist.Hist(
                  "Events",
                  hist.Cat("sample","XCone Jet p_t"),
                  hist.Bin("jet_pt","XCone Jet p_t", 90, 0, 1000)),
            "jet_eta": hist.Hist(
                  "Events",
                  hist.Cat("sample","XCone Jet eta"),
                  hist.Bin("jet_eta","XCone Jet eta", 90, -3, 3)),
            "jet_phi": hist.Hist(
                  "Events",
                  hist.Cat("sample","XCone Jet phi"),
                  hist.Bin("jet_phi","XCone Jet phi", 90, -3.5, 3.5)),
            "jet_high_pt": hist.Hist(
                  "Events",
                  hist.Cat("sample","XCone Jet high pt"),
                  hist.Bin("jet_high_pt","XCone Jet high pt", 90, 0, 1000)),
        })


    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    def process(self, events):
        dataset = events.metadata['dataset'] 

        particles = events.GenPart.pdgId
        values = ak.flatten(particles)
		

        e = events.Electron
        j = events.XConeJet
        j = ak.with_name(j, "PtEtaPhiMLorentzVector")
        mu = events.Muon
        print("Total events: "+ str(ak.num(e, axis=0)))
        
        #Electron tight cut
        e['isPresElec'] = isPresElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all)
        e['isLooseElec'] = isLooseElec(e.miniPFRelIso_all, e.sip3d, e.lostHits)
        e = e[e.isPresElec & e.isLooseElec]
        elecpt = ak.flatten(e.pt)
        elecCount = ak.num(e, axis=1)
        print("Total events after electron cut: " + str(ak.count_nonzero(elecCount)))
        
        #Muon tight cut
        mu['isPresMuon']= isPresMuon(mu.dxy, mu.dz, mu.eta, mu.pt)
        mu['isLooseMuon'] = isLooseMuon(mu.miniPFRelIso_all, mu.sip3d)
        mu = mu[mu.isPresMuon & mu.isLooseMuon]
        muonpt = ak.flatten(mu.pt)
        muonCount = ak.num(mu,axis=1)
        print("Total events after muon cut: " + str(ak.count_nonzero(muonCount)))
        
        #Single Tight Lepton
        singLep = elecCount + muonCount ==1
        print("Total events after electron and muon cut: " + str(ak.count_nonzero(singLep)))

        j['isClean'] = isClean(j,e,drmin=0.4) & isClean(j,mu,drmin=0.4)
        j = j[j.isClean]
        j = j[singLep]
        j = j[j.pt>=25]
        j = j[abs(j.eta)<2.5]
        jetCount = ak.num(j,axis=1)
        j = j[jetCount>=4]
        print("Total events after single lepton and jetpt cut: " + str(ak.count_nonzero(jetCount>=4)))
        
        #Two loose btags, one medium.
        btagLoose = j.btagDeepFlavB
        btagLoose = btagLoose[btagLoose>0.0490]
        btagCountLoose = ak.num(btagLoose, axis = 1)
        j = j[btagCountLoose>=2]
        btagMed = j.btagDeepFlavB
        btagMed = btagMed[btagMed>0.2783]
        btagCountMed = ak.num(btagMed, axis=1)
        j = j[btagCountMed>=1]
        print("Total events after lep, jetpt, and btag cut: " + str(ak.count_nonzero(btagCountMed>=1)))
        jetpt = ak.flatten(j.pt)
        jeteta = ak.flatten(j.eta)
        jetphi = ak.flatten(j.phi)
        jethighpt = ak.max(j.pt, axis=1)

        # fill Histos
        hout = self.accumulator.identity()
        hout['events'].fill(
            sample=dataset,
            pdgID=values,
        )
        hout['electron_pt'].fill(
            sample = dataset,
            electron_pt = elecpt,
        )
        hout['muon_pt'].fill(
            sample = dataset,
            muon_pt = muonpt,
        )
        hout['jet_pt'].fill(
            sample = dataset,
            jet_pt = jetpt,
        )
        hout['jet_eta'].fill(
            sample = dataset,
            jet_eta = jeteta,
        )
        hout['jet_phi'].fill(
            sample = dataset,
            jet_phi = jetphi,
        )
        hout['jet_high_pt'].fill(
        	sample = dataset,
        	jet_high_pt = jethighpt,
        )
        return hout

    def postprocess(self, accumulator):
        return accumulator


