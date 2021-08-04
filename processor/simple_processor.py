#!/usr/bin/env python
import coffea
import numpy as np
import awkward as ak
from coffea import hist, processor
#from topcoffea.modules.objects import isTightElec, isTightMuon
#from coffea.analysis_tools import PackedSelection

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
		
		#Electron tight cut
        e = events.Electron
        print("Total events: "+ str(ak.num(e, axis=0)))
        
        e['isTightElec'] = isTightElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all, e.sip3d, e.mvaTTH, e.mvaFall17V2Iso, e.lostHits, e.convVeto, e.tightCharge, e.sieie, e.hoe, e.eInvMinusPInv, minpt=30)
        e = e[e.isTightElec]
        elec = ak.flatten(e.pt)
        elecCount = ak.num(e, axis=1)
        print("Total events after electron cut: " + str(ak.count_nonzero(elecCount)))
        
        #Muon tight cut
        mu = events.Muon
        mu['isTightMuon']= isTightMuon(mu.pt, mu.eta, mu.dxy, mu.dz, mu.pfRelIso03_all, mu.sip3d, mu.mvaTTH, mu.mediumPromptId, mu.tightCharge, mu.looseId, minpt=30)
        mu = mu[mu.isTightMuon]
        muonCount = ak.num(mu,axis=1)
        print("Total events after muon cut: " + str(ak.count_nonzero(muonCount)))
        
        #Single Tight Lepton
        singLep = elecCount + muonCount ==1
        print("Total events after electron and muon cut: " + str(ak.count_nonzero(singLep)))

        j = events.XConeJet.pt
        #j = j[(abs(e.pt)>=30) & (abs(e.eta)<=2.5)]
        j = j[singLep]
        jet = ak.flatten(j)
        #print("Total events after electron and muon cut: " + str(ak.num(ak.flatten(singLep), axis=0)))
            

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


