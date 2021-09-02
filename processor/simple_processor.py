#!/usr/bin/env python
import coffea
import numpy as np
import awkward as ak
from coffea import hist, processor
import topcoffea.modules.objects as top


def isTightElec(pt, eta, dxy, dz, miniIso, sip3D, mvaTTH, elecMVA, lostHits, convVeto, tightCharge, sieie, hoe, eInvMinusPInv, minpt=15.0):
  maskPOGMVA = ((pt<10)&(abs(eta)<0.8)&(elecMVA>-0.13))|((pt<10)&(abs(eta)>0.8)&(abs(eta)<1.44)&(elecMVA>-0.32))|((pt<10)&(abs(eta)>1.44)&(elecMVA>-0.08))|\
               ((pt>10)&(abs(eta)<0.8)&(elecMVA>-0.86))|((pt>10)&(abs(eta)>0.8)&(abs(eta)<1.44)&(elecMVA>-0.81))|((pt>10)&(abs(eta)>1.44)&(elecMVA>-0.72))
  maskSieie  = ((abs(eta)<1.479)&(sieie<0.011))|((abs(eta)>1.479)&(sieie<0.030))
  maskhoe    = ((abs(eta)<1.479)&(hoe<0.10))|((abs(eta)>1.479)&(hoe<0.07))
  mask = (pt>minpt)&(abs(eta)<2.5)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(sip3D<8)&(lostHits<=1)&\
         (convVeto)&(maskSieie)&(maskPOGMVA)&(eInvMinusPInv>-0.04)&(maskhoe)&(miniIso<0.25)&(mvaTTH>0.90)&(tightCharge==2)
  return mask
  
def isTightMuon(pt, eta, dxy, dz, miniIso, sip3D, mvaTTH, mediumPrompt, tightCharge, looseId, minpt=10.0):
  mask = (pt>minpt)&(abs(eta)<2.5)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(sip3D<8)&(looseId)&(miniIso<0.25)&(mvaTTH>0.90)&(tightCharge==2)&(mediumPrompt)
  return mask

def deltaR(eta1, phi1, eta2, phi2):
	return np.sqrt((eta1-eta2)**2 + (np.arccos(np.cos(phi1))-np.arccos(np.cos(phi2)))**2)
	

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self):

        # Create the histograms
        # In general, histograms depend on 'sample', 'channel' (final state) and 'cut' (level of selection)
        self._accumulator = processor.dict_accumulator({
            "events": hist.Hist(
                  "Events",
                  hist.Cat("sample", "A"),
                  hist.Bin("pdgID", "B", 1200, -600, 600 )),
            "electrons": hist.Hist(
                  "Electrons",
                  hist.Cat("sample", "A"),
                  hist.Bin("pt", "B", 600, 0, 600)),
            "jets": hist.Hist(
                  "Jets P_T",
                  hist.Cat("sample","A"),
                  hist.Bin("jets","Events", 90, 0, 1000)),
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
        mu = events.Muon
        print("Total events: "+ str(ak.num(e, axis=0)))
        
        #Electron tight cut
        e['isTightElec'] = isTightElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all, e.sip3d, e.mvaTTH, e.mvaFall17V2Iso, e.lostHits, e.convVeto, e.tightCharge, e.sieie, e.hoe, e.eInvMinusPInv, minpt=30)
        e = e[e.isTightElec]
        #e['isPres'] = top.isPresElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all, e.sip3d, getattr(e,"mvaFall17V2noIso_WPL"))
        #e['isLooseE'] = top.isLooseElec(e.miniPFRelIso_all,e.sip3d,e.lostHits)
        #e = e[e.isPres & e.isLooseE]
        elec = ak.flatten(e.pt)
        elecCount = ak.num(e, axis=1)
        print("Total events after electron cut: " + str(ak.count_nonzero(elecCount)))
        
        #Muon tight cut
        mu['isTightMuon']= isTightMuon(mu.pt, mu.eta, mu.dxy, mu.dz, mu.pfRelIso03_all, mu.sip3d, mu.mvaTTH, mu.mediumPromptId, mu.tightCharge, mu.looseId, minpt=30)
        mu = mu[mu.isTightMuon]
        #mu['isPres'] = top.isPresMuon(mu.dxy, mu.dz, mu.sip3d, mu.eta, mu.pt, mu.miniPFRelIso_all)
        #mu['isLooseM'] = top.isLooseMuon(mu.miniPFRelIso_all,mu.sip3d,mu.looseId)
        #mu = mu[mu.isPres & mu.isLooseM]
        muonCount = ak.num(mu,axis=1)
        print("Total events after muon cut: " + str(ak.count_nonzero(muonCount)))
        
        #Single Tight Lepton
        singLep = elecCount + muonCount ==1
        print("Total events after electron and muon cut: " + str(ak.count_nonzero(singLep)))


        j = j[singLep]
        j = j[j.pt>=25]
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
        jet = ak.flatten(j.pt)

        # fill Histos
        hout = self.accumulator.identity()
        hout['events'].fill(
            sample=dataset,
            pdgID=values,
        )
        hout['electrons'].fill(
            sample = dataset,
            pt = elec,
        )
        hout['jets'].fill(
            sample = dataset,
            jets = jet,
        )
        return hout

    def postprocess(self, accumulator):
        return accumulator


