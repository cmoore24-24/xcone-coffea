#!/usr/bin/env python
import coffea
import numpy as np
import awkward as ak
from coffea import hist, processor
import topcoffea.modules.objects as top

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
        e['isPres'] = top.isPresElec(e.pt, e.eta, e.dxy, e.dz, e.miniPFRelIso_all, e.sip3d, getattr(e,"mvaFall17V2noIso_WPL"))
        e['isLooseE'] = top.isLooseElec(e.miniPFRelIso_all,e.sip3d,e.lostHits)
        e = e[e.isPres & e.isLooseE]
        elec = ak.flatten(e.pt)
        elecCount = ak.num(e, axis=1)
        print("electron cut: " + str(ak.count_nonzero(elecCount)))
 
        #Muon tight cut
        mu = events.Muon
        mu['isPres'] = top.isPresMuon(mu.dxy, mu.dz, mu.sip3d, mu.eta, mu.pt, mu.miniPFRelIso_all)
        mu['isLooseM'] = top.isLooseMuon(mu.miniPFRelIso_all,mu.sip3d,mu.looseId)
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


