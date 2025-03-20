import awkward as ak
from coffea import processor
from coffea.nanoevents.methods import candidate
#import hist
from coffea import hist
import numpy as np
from HMTprocessor import util
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods import vector
ak.behavior.update(candidate.behavior)

from HMTprocessor.BaseMDSnanoProcessor import BaseMDSnanoProcessor

class MDStrigProcessor(BaseMDSnanoProcessor):
    """
    Processor for MDS nano for  muon shower in the muon system. 

    Parameters 
    ----------
    """

    @property
    def accumulator(self):
        return self._accumulator

    def __init__(self,**options):
        super().__init__()
        defaultOptions = { 'forTrig':False}
        options = { **defaultOptions, **options }

        self._forTrig = options['forTrig']
        histograms={}
        histograms['sumw']= processor.defaultdict_accumulator(float)
        histograms['noiseFilters']= processor.defaultdict_accumulator(float)
        self._accumulator = processor.dict_accumulator( histograms )
        self.isData = True
        pass    

    
    
    def process(self, events):
        
        dataset = events.metadata['dataset']

        ## apply golden json
        events = self.ApplyGoldenJSON(events)
        ## add noise filters
        events.PassNoiseFilters = self.AddNoiseFilters(events)
       
        clsName = "cscMDSHLTCluster" if "cscMDSHLTCluster" in events.fields else "MDSHLTCluster"
        cls = getattr(events,clsName)
        cls = self.build_cls(cls)
        muons = events.Muon
        jets = events.Jet

        hlt   = events.HLT        
        L1   = events.L1
        output = self.accumulator.identity() 

        hall = self.fillEff(hlt,L1,muons,jets,cls,events.PassNoiseFilters,dataset)
        #print(hall)
        for key,h in hall.items():
            output[key] = h        
        hall = self.fill_rateMap(cls,events)
        for key,h in hall.items():
            output[key] = h        

        ## output hist
        output["sumw"][dataset] += len(events)
        output["noiseFilters"][dataset] += ak.sum(events.PassNoiseFilters)

        return output

    def fill_rateMap(self,cls,events):

        hall = {}
        runMin = 378971
        runMax = 386951 

        hall["Reg_runNumber"] = hist.Hist("Events",hist.Cat("reg","reg"),hist.Bin("runNumber", "runNumber", 100,runMin,runMax))     

        hall["cls_etaPhi"] = hist.Hist("Events",hist.Bin("eta", "eta", 30,-3,3), hist.Bin("phi", "phi", 32,0,3.2))  
        hall["cls_RZ"] = hist.Hist("Events",hist.Bin("R", "R", 70,0,700), hist.Bin("z", "z", 100,-1100,1100))  
        hall["nCluster"] = hist.Hist("Events",hist.Bin("nCls", "nCls", 10,0,10))  

        evt_sel = events.L1.SingleMuShower_Nominal_pRECO 
        for sign in ["+","-"]:
            for regName,region in self.getRegions(cls).items():
                if sign =="+":
                    cls_sel = ak.any(cls[region].z>0,axis=1)
                else:
                    cls_sel = ak.any(cls[region].z<0,axis=1)            
                hall["Reg_runNumber"].fill(reg=sign+regName,runNumber = events.run[evt_sel & cls_sel])

        hall["cls_etaPhi"].fill(eta=ak.flatten(cls[evt_sel].eta),
                                phi=ak.flatten(cls[evt_sel].phi))
        hall["cls_RZ"].fill(R=ak.flatten(cls[evt_sel].R),
                            z=ak.flatten(cls[evt_sel].z))
        hall["nCluster"].fill(nCls=ak.num(cls[evt_sel]))
        
        return hall
 
    def fillEff(self,hlt,L1,muons,jets,cls,noise,dataset):
        hall={}
         
        #good muons
        good_muons = muons[(muons.pfRelIso04_all<0.4)&(muons.pt>20)&(abs(muons.eta)<2.4)&(muons.looseId)]
        tight_muons = muons[(muons.pfRelIso04_all<0.4)&(muons.pt>20)&(abs(muons.eta)<2.4)&(muons.tightId)]
        # cluster matched to any good muons
        cls  = util.match(cls,good_muons,"matched_muon",0.4)
        cls  = util.match(cls,tight_muons,"matched_Tightmuon",0.4)
        #lepton jets
        lepJets = jets[(jets.jetId&3)!=3]

        # cluster matched to any bad muons
        cls  = util.match(cls,muons,"matched_anymuon",0.4)
        cls  = util.match(cls,jets,"matched_jets",0.4)
        cls  = util.match(cls,lepJets,"matched_lepjets",0.4)

        cls_sel           = ((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(cls.nStation==1)&(cls.matched_muon==1)
        cls_sel_noLepJet  = ((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(cls.nStation==1)&(cls.matched_muon==1)&(cls.matched_lepjets==0)
        cls_sel_tightMu   = ((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(cls.nStation==1)&(cls.matched_Tightmuon==1)
        cls_sel_noMu      = ((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(cls.nStation==1)&(cls.matched_anymuon==0)
        cls_sel_noMuNoJet = ((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(cls.nStation==1)&(cls.matched_anymuon==0)&(cls.matched_jets==0)

        #cls_sel_nChamber1  = ((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(cls.nStation==1)&(cls.matched_muon==1)&(cls.nChambers==1)
        #cls_sel_nChamberGe1  = ((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(cls.nStation==1)&(cls.matched_muon==1)&(cls.nChambers>=1)

        #Fill diagnosis plots
        regMap =  self.getRegions(cls)
        evt_sel = (ak.num(cls)==1)
        evt_sel_noNoise = (ak.num(cls)==1) & (noise)

        ## Hotspot filter
        cls_noiseVeto = ~( (cls.z<0)&(regMap["ME31"] | regMap["ME41"]) & (((cls.phi>0.4) & (cls.phi<0.8))|(abs(cls.phi)>2.8)))
        cls_sel  = (cls_sel) & (cls_noiseVeto)
        
        #cluster prop of denom     
        denom = (cls_sel) & (evt_sel) & (hlt.IsoMu24==1) 
        hlist = self.fillbasic({(f"IsoMu24_denom","all"):cls[denom]},None)
        denom = (cls_sel_noMuNoJet) & (evt_sel) & (hlt.IsoMu24==1)
        hlist = self.fillbasic({(f"IsoMu24_unMatch_denom","all"):cls[denom]},hlist)
        
        # Large cluster properties
        for reg in ["ME21","ME31","ME41"]:
            denom = (cls_sel)&(cls.size>=400) & (evt_sel) & (regMap[reg]) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            hlist = self.fillbasic(
                        {(f"TrgMatchMu_numer",reg):cls[numer],(f"TrgMatchMu_denom",reg):cls[denom]},hlist
            )
            denom = (cls_sel)&(cls.size>=400) & (evt_sel_noNoise) & (regMap[reg]) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            ## Fill no noise
            hlist = self.fillbasic(
                        {(f"TrgMatchMuNoNoise_numer",reg):cls[numer],(f"TrgMatchMuNoNoise_denom",reg):cls[denom]},
                        hlist
            )
            ## Fill unMatched 
            denom = (cls_sel_noMuNoJet)&(cls.size>=400) & (evt_sel_noNoise) & (regMap[reg]) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            hlist = self.fillbasic(
                        {(f"TrgUnMatch_numer",reg):cls[numer],(f"TrgUnMatch_denom",reg):cls[denom]},hlist
            )
       
        for hname,h in hlist.items():
            key = "_".join([ ax.name for ax in h.dense_axes()])
            hall[key] = h        

        bin_size = [   0.,  50.,   60., 70,  80.,90,  100.,110,  120.,130,  140.,150,  160.,
                     180.,190, 200.,  220.,  240.,  260.,  280.,  300.,  320.,  340.,
                     360., 380.,  400.,  420.,  440.,  460.,  480.,  500.,  520.,
                     540., 560.,  580.,  600.,  620.,  640.,  660.,  680.,  700.,
                     720., 740.,  760.,  780.,  800.,  820.,  840.,  860.,  880.,
                     900., 920.,  940.,  960.,  980., 1000.,  1500,  2000, 4000.  ]

        # book histograms
        for hname in [
            "eff_IsoMu24_Cls100Mu_reg",
            "eff_Ele30_Cls100Ele_reg",
            "eff_IsoMu24_ClsLoose",
            "eff_IsoMu24_ClsLoose_L1tight",
            "eff_IsoMu24_ClsLoose_noPhi",
            "eff_IsoMu24_ClsLoose_noLepJet",
            "eff_IsoMu24_ClsLoose_noNoise",
            "eff_IsoMu24_ClsLoose_tightMu",
            "eff_IsoMu24_ClsLoose_noMatchMu",
            "eff_IsoMu24_ClsLoose_unMatch",
            "eff_IsoMu24_ClsLoose_pos",
            "eff_IsoMu24_ClsLoose_neg",
            #"eff_IsoMu24_ClsLoose_nCham1",
            #"eff_IsoMu24_ClsLoose_nChamGe1"
            ]:
            hall[hname] = hist.Hist("Events",hist.Cat("dataset","dataset"),
                                           hist.Cat("sample","sample"),hist.Cat("reg","reg"),
                                           hist.Bin("ClusterSize", "ClusterSize", bin_size)) 

        # Fill matchLepJet
        hall["h_matchedLepJet"] = hist.Hist("Events",hist.Cat("dataset","dataset"),
                                           hist.Cat("sample","sample"),hist.Cat("reg","reg"),
                                           hist.Bin("matchedLepJet", "matchedLepJet", 2,0,2))  
        
        for regName,region in self.getRegions(cls).items():
            denom = (cls_sel) & (evt_sel) & (region) & (hlt.IsoMu24==1)       
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            hall["h_matchedLepJet"].fill(sample="denom",matchedLepJet=self.prep(cls[denom].matched_lepjets),reg=regName,dataset=dataset)
            hall["h_matchedLepJet"].fill(sample="numer",matchedLepJet=self.prep(cls[numer].matched_lepjets),reg=regName,dataset=dataset)
 
        ## Fill nChambers
        #hall["h_nChamber"] = hist.Hist("Events",hist.Cat("dataset","dataset"),
        #                                   hist.Cat("sample","sample"),hist.Cat("reg","reg"),
        #                                   hist.Bin("nChambers", "nChambers", 40,0,40))  
        #
        #for regName,region in self.getRegions(cls).items():
        #    denom = (cls_sel) & (evt_sel) & (region) & (hlt.IsoMu24==1)       
        #    hall["h_nChamber"].fill(sample="matchedMuon",nChambers=self.prep(cls[denom].nChambers),reg=regName,dataset=dataset)
        #    denom = (cls_sel_noMu) & (evt_sel) & (region) & (hlt.IsoMu24==1)       
        #    hall["h_nChamber"].fill(sample="NoMuon",nChambers=self.prep(cls[denom].nChambers),reg=regName,dataset=dataset)
        #    denom = (cls_sel_noMuNoJet) & (evt_sel) & (region) & (hlt.IsoMu24==1)       
        #    hall["h_nChamber"].fill(sample="Unmatched",nChambers=self.prep(cls[denom].nChambers),reg=regName,dataset=dataset)

        # Fill Efficiencies {"ME11":cls within ME11}
        for regName,region in self.getRegions(cls).items():
            denom = (cls_sel) & (evt_sel) & (region) & (hlt.IsoMu24==1)       
            numer = (denom)&(hlt.CscCluster100_Mu5==1)

            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_Cls100Mu_reg"],
                       reg=regName,dataset=dataset)
            
            denom = (cls_sel) & (evt_sel) & (region) & (hlt.Ele30_WPTight_Gsf==1)
            numer = (denom)&(hlt.CscCluster100_Ele5==1)

            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_Ele30_Cls100Ele_reg"],
                       reg=regName,dataset=dataset)

            denom = (cls_sel) & (evt_sel) & (region) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose"],
                       reg=regName,
                       dataset=dataset
                    )
            numer = (denom)&(L1.SingleMuShower_Tight_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_L1tight"],
                       reg=regName,
                       dataset=dataset
                    )

            denom = (cls_sel_noMu) & (evt_sel) & (region) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_noMatchMu"],
                       reg=regName,dataset=dataset
                    )
            denom = (cls_sel_noMuNoJet) & (evt_sel) & (region) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_unMatch"],
                       reg=regName,dataset=dataset
                    )
            denom = (cls_sel) & (evt_sel) & (region) & ~((cls.phi<0.6)&(cls.phi>0.0)) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_noPhi"],
                       reg=regName,dataset=dataset
                    )
            denom = (cls_sel) & (evt_sel) & (region) & (cls.z>0) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_pos"],
                       reg=regName,dataset=dataset
                    )
            denom = (cls_sel) & (evt_sel) & (region) & (cls.z<0) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_neg"],
                       reg=regName,dataset=dataset
                    )
            denom = (cls_sel_tightMu) & (evt_sel) & (region) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_tightMu"],
                       reg=regName,dataset=dataset
                    )
            denom = (cls_sel) & (evt_sel_noNoise) & (region) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_noNoise"],
                       reg=regName,dataset=dataset
                    )
            denom = (cls_sel_noLepJet) & (evt_sel) & (region) & (hlt.IsoMu24==1)
            numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            self.addEff(self.prep(cls[numer].size),
                        self.prep(cls[denom].size),
                       hall["eff_IsoMu24_ClsLoose_noLepJet"],
                       reg=regName,dataset=dataset
                    )

            #denom = (cls_sel_nChamber1) & (evt_sel) & (region) & (cls.z<0) & (hlt.IsoMu24==1)
            #numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            #self.addEff(self.prep(cls[numer].size),
            #            self.prep(cls[denom].size),
            #           hall["eff_IsoMu24_ClsLoose_nCham1"],
            #           reg=regName,
            #           dataset=dataset
            #        )
            #denom = (cls_sel_nChamberGe1) & (evt_sel) & (region) & (cls.z<0) & (hlt.IsoMu24==1)
            #numer = (denom)&(L1.SingleMuShower_Nominal_pRECO==1)
            #self.addEff(self.prep(cls[numer].size),
            #            self.prep(cls[denom].size),
            #           hall["eff_IsoMu24_ClsLoose_nChamGe1"],
            #           reg=regName,
            #           dataset=dataset
            #        )

        return hall  

    def addEff(self,numer,denom,h=None,**cats):
        if h==None:
            h= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("ClusterSize", "ClusterSize", 50, 0, 1000))        
        h.fill(sample="numer",ClusterSize=numer,**cats)
        h.fill(sample="denom",ClusterSize=denom,**cats)
        return h    
    
    def postprocess(self, accumulator):
        pass
        
