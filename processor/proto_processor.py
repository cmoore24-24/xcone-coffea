#!/usr/bin/env python
import coffea
import numpy as np
import awkward as ak
from coffea import hist, processor
import topcoffea.modules.objects as top

def isPresElec(pt, eta, dxy, dz, miniIso):
    mask = (pt>32)&(abs(eta)<2.4)&(abs(dxy)<0.05)&(abs(dz)<0.1)&(miniIso<0.4)
    return mask
def isLooseElec(miniPFRelIso_all,sip3d,lostHits):
    return (miniPFRelIso_all<0.4) & (sip3d<8) & (lostHits<=1)
def isPresMuon(dxy, dz, eta, pt):
    mask = (abs(dxy)<0.2)&(abs(dz)<0.5)&(abs(eta)<2.4)&(pt>27)
    return mask
def isLooseMuon(miniPFRelIso_all,sip3d):
    return (miniPFRelIso_all<0.4) & (sip3d<8)


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
        print("Total : " + str(ak.num(e, axis=0)))
        e['isPres'] = isPresElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all)
        e['isLooseE'] = isLooseElec(e.miniPFRelIso_all,e.sip3d,e.lostHits)
        e = e[e.isPres & e.isLooseE]
        elec = ak.flatten(e.pt)
        elecCount = ak.num(e, axis=1)
        print("electron cut: " + str(ak.count_nonzero(elecCount)))
 
        #Muon tight cut
        mu = events.Muon
        mu['isPres'] = isPresMuon(mu.dxy, mu.dz, mu.eta, mu.pt)
        mu['isLooseM'] = isLooseMuon(mu.miniPFRelIso_all,mu.sip3d)
        mu = mu[mu.isPres & mu.isLooseM]
        muonCount = ak.num(mu,axis=1)
        print("muon cut: " + str(ak.count_nonzero(muonCount)))





        #Single Tight Lepton
        singLep = elecCount + muonCount ==1
        print("total count: " + str(ak.count_nonzero(singLep)))

        j = events.XConeJet.pt
        #j = j[(abs(e.pt)>=30) & (abs(e.eta)<=2.5)]
        j = j[singLep]
        jet = ak.flatten(j)
        
            

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

# Electron
# pt > 35
# |eta| < 1.442 OR 1.566 < |eta| < 2.4 ( |eta| < 2.4  ?)




#Muon
# Pt > 29
# |eta| < 2.4
# |dxy| < 0.2
# |dz| < 0.5
# sip3d < 4,8

