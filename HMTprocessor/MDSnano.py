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

class MDSnanoProcessor(processor.ProcessorABC):
    """
    Processor for MDS nano for  muon shower in the muon system. 

    Parameters 
    ----------
    forTrig : bool (default is False)
        Run trigger histograms 
    """

    def __init__(self,**options):
        defaultOptions = { 'forTrig':False}
        options = { **defaultOptions, **options }
    
        self._forTrig = options['forTrig']
        histograms={}
        histograms['sumw']= processor.defaultdict_accumulator(float)
        self._accumulator = processor.dict_accumulator( histograms )
        pass

    def fillbasic(self,samples={}):
        h1= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("nClusters", "nClusters", 10, 0, 10))
        h2= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("eta", "eta", 40, -5, 5))    
        h3= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("phi", "phi", 40, -np.pi, np.pi))        
        h4= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("size", "size", 100, 0, 1000))        
        h5= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("time", 'ClusterTime[ns]',80, -80, 80))            
        h6 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("nME11_12", 'nME11_12',40, 0, 40))        
        h12= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("timeSpread", "timeSpread[ns]", 80, 0, 100))
            
        h7 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("x", "x", 50, -1000, 1000))        
        h8 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("y", "y", 50, -1000, 1000))        
        h9 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("z", "z", 50, -1000, 1000))            
        h13 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("R", "R", 50, 0, 1000))            
        h10= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("Nstation10", "Nstation10", 8, 0, 8))    
        h11= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("AvgStation10", "AvgStation10", 40, 0, 8))  

        hall = [h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13]
        for dataset,cls in samples.items():
            h1.fill(sample=dataset,nClusters=ak.num(cls))
            h2.fill(sample=dataset,eta = ak.firsts(cls.eta))
            h3.fill(sample=dataset,phi = ak.firsts(cls.phi)) 
            h4.fill(sample=dataset,size = ak.firsts(cls.size))
            h5.fill(sample=dataset,time = ak.firsts(cls.time)) 
            h6.fill (sample=dataset,nME11_12 = ak.firsts(cls.nME11+cls.nME12))
            h12.fill(sample=dataset,timeSpread = ak.firsts(cls.timeSpread))                
        
            h7.fill( sample=dataset,x = ak.firsts(cls.x))    
            h8.fill( sample=dataset,y = ak.firsts(cls.y))    
            h9.fill( sample=dataset,z = ak.firsts(cls.z))
            h13.fill( sample=dataset,R = ak.firsts((cls.x**2+cls.y**2)**0.5))
            h10.fill(sample=dataset,Nstation10 = ak.firsts(cls.nStation))            
            h11.fill(sample=dataset,AvgStation10 = ak.firsts(cls.avgStation)) 

        return hall 

    def getRegions(self,cls):    
        ME11 = (cls.R>100)&(cls.R<275) &(abs(cls.z)>580)&(abs(cls.z)<632) 
        ME12 = (cls.R>275)&(cls.R<465) &(abs(cls.z)>668)&(abs(cls.z)<724)
        ME13 = (cls.R>505)&(cls.R<700) &(abs(cls.z)>668)&(abs(cls.z)<724)
    
        ME21 = (cls.R>139)&(cls.R<345) &(abs(cls.z)>789)&(abs(cls.z)<850)
        ME22 = (cls.R>357)&(cls.R<700) &(abs(cls.z)>791)&(abs(cls.z)<850)
    
        ME31 = (cls.R>160)&(cls.R<345) &(abs(cls.z)>915)&(abs(cls.z)<970)
        ME32 = (cls.R>357)&(cls.R<700) &(abs(cls.z)>911)&(abs(cls.z)<970)
    
        ME41 = (cls.R>178)&(cls.R<345) &(abs(cls.z)>1002)&(abs(cls.z)<1063)
        ME42 = (cls.R>357)&(cls.R<700) &(abs(cls.z)>1002)&(abs(cls.z)<1063)
    
        regions = {        "ME11":ME11,         "ME12":ME12,         "ME13":ME13, 
                            "ME21":ME21,         "ME22":ME22, 
                            "ME31":ME31,         "ME32":ME32, 
                            "ME41":ME41,         "ME42":ME42, 
        }
        return regions

    def addEff(self,numer,denom,h=None,**cats):
        if h==None:
            h= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("ClusterSize", "ClusterSize", 50, 0, 1000))        
        h.fill(sample="numer",ClusterSize=numer,**cats)
        h.fill(sample="denom",ClusterSize=denom,**cats)
        return h


    def fillEff(self,hlt,cls):
        hall={}
        #denom = hlt.Mu50
        #denom = (hlt.Mu50)&(ak.any(((cls.nME11+cls.nME12)==0),axis=1))        
        #denom = (hlt.Mu50)&(ak.any(((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5),axis=1)
         
        cls = ak.firsts(cls)
        denom = (hlt.IsoMu24)&((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)         
        numer = (denom)&(hlt.CscCluster100_Mu5==1)        

        hall["eff_IsoMu24_Cls100Mu"] = self.addEff(ak.fill_none(cls[numer].size,-1),
                                                   ak.fill_none(cls[denom].size,-1))

        hall["eff_IsoMu24_Cls100Mu_reg"] = hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),
                                           hist.Bin("ClusterSize", "ClusterSize", 50, 10, 1000))        
        
        hall["eff_Ele30_Cls100Ele_reg"] = hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),
                                           hist.Bin("ClusterSize", "ClusterSize", 50, 10, 1000))        
        hall["eff_IsoMu24_ClsLoose"] = hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("cat","cat"),hist.Cat("reg","reg"),
                                           hist.Bin("ClusterSize", "ClusterSize", 50, 10, 1000))        
        for regName,region in self.getRegions(cls).items():
            denom = (hlt.IsoMu24)&((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(region)
            numer = (denom)&(hlt.CscCluster100_Mu5==1)

            self.addEff(ak.fill_none(cls[numer].size,-1),
                        ak.fill_none(cls[denom].size,-1),
                       hall["eff_IsoMu24_Cls100Mu_reg"],
                       reg=regName)
            
            denom = (hlt.Ele30_WPTight_Gsf)&((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(region)
            numer = (denom)&(hlt.CscCluster100_Ele5==1)

            self.addEff(ak.fill_none(cls[numer].size,-1),
                        ak.fill_none(cls[denom].size,-1),
                       hall["eff_Ele30_Cls100Ele_reg"],
                       reg=regName)

            denom = (hlt.IsoMu24)&((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(region)&(cls.nStation==1)
            numer = (denom)&(hlt.CscCluster_Loose==1)
            self.addEff(ak.fill_none(cls[numer].size,-1),
                        ak.fill_none(cls[denom].size,-1),
                       hall["eff_IsoMu24_ClsLoose"],
                       reg=regName,cat = "nStat_eq1"
                    )
            denom = (hlt.IsoMu24)&((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5)&(region)&(cls.nStation>1)
            numer = (denom)&(hlt.CscCluster_Loose==1)
            self.addEff(ak.fill_none(cls[numer].size,-1),
                        ak.fill_none(cls[denom].size,-1),
                       hall["eff_IsoMu24_ClsLoose"],
                       reg=regName,cat = "nStat_gt1"
                    )

            
        return hall
        
    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, events):

        dataset = events.metadata['dataset']

        cls   = events.MDSHLTCluster
        cls = ak.with_field(cls,(cls.x**2+cls.y**2)**0.5,"R")
        ## sort cls in event according to nHit
        idx   = ak.argsort(cls.size,ascending=False)
        cls   = cls[idx] 
        muons = events.Muon
        hlt   = events.HLT        
        binning = np.array(
            [ 0.,  20.,  40.,  60.,  80., 100., 120., 140., 160., 180., 200.,
           220., 240., 260., 280., 300., 500., 700.]
        )

        output = self.accumulator.identity() 
      
        if self._forTrig:
            denom = (hlt.IsoMu24)&(ak.any(((cls.nME11+cls.nME12)==0)&(cls.time<12.5)&(cls.time>-5),axis=1))
            numer = (denom)&(hlt.CscCluster100_Mu5==1)
    
            hall = self.fillbasic(
                        {"numer":cls[numer],
                         "denom":cls[denom]
                        }
            )
            for h in hall:
                key = "Trg_"+h.dense_axes()[0].name
                output[key] = h        
    
            hall = self.fillEff(hlt,cls)
            print(hall)
            for key,h in hall.items():
                output[key] = h
        else:
            samples={
                "Trigger":events.MDSHLTCluster[hlt.CscCluster100_Mu5]   
            }
            hall = fillbasic(samples)
            h3= hist.Hist("Events",hist.Bin("run","run",100,378981,380470),hist.Bin("phi", "phi", 32, -np.pi, np.pi))
            h3.fill(run=events.run,phi= ak.fill_none(ak.firsts(cls.phi),-999))
            output["h_run_phi"] = h3

            h3= hist.Hist("Events",hist.Bin("eta", "eta", 60, -3, 3),hist.Bin("phi", "phi", 64, -np.pi, np.pi)) 
            h3.fill(eta=ak.fill_none(ak.firsts(cls.eta),-999),phi= ak.fill_none(ak.firsts(cls.phi),-999))
            output['h_eta_phi'] = h3

        ## output hist
        output["sumw"][dataset] += len(events)
         
        return output


    def postprocess(self, accumulator):
        pass
