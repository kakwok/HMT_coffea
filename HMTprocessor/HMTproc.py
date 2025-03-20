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

class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        histograms={}
        histograms['sumw']= processor.defaultdict_accumulator(float)
        self._accumulator = processor.dict_accumulator( histograms )
        pass

    @property
    def accumulator(self):
        return self._accumulator

    def pack(self,events):
        events.info=ak.zip({
            "runNum": events.runNum,
            "lumiNum": events.lumiNum,
            "evtNum": events.eventNum,
            }
        )
        events.cls = ak.zip(
            {k.replace("ca4CSCcluster",""):getattr(events,k) for k in events.fields if k.startswith("ca4CSCcluster")}
            ,with_name="PtEtaPhiMLorentzVector", 
            behavior=vector.behavior
            )
        
        #events.muons = ak.zip(
        #    {k.replace("muon",""):getattr(events,k) for k in events.fields if k.startswith("muon")}
        #    ,with_name="PtEtaPhiMLorentzVector", 
        #    behavior=vector.behavior
        #    )
        
        return events

    def fillHMT(self,lctHMT,output,tag):

        #labels  = ["ME11","ME12","ME13",'ME21','ME22','ME31','ME32',"ME41","ME42"]
        #sr_map  = [(8,9),(7,10),(6,11),(5,12),(4,13),(3,14),(2,15),(1,16),(0,17)]

        labels  = ["ME13",'ME21','ME22','ME31','ME32',"ME41","ME42"]
        sr_map  = [(6,11),(5,12),(4,13),(3,14),(2,15),(1,16),(0,17)]
        for i,label in enumerate(labels):
            #h = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("size", "size", 20, 0, 200))
            #pick station-ring
            #sel = ((lctHMT.sr==sr_map[i][0])|lctHMT.sr==sr_map[i][1])
            #h.fill(sample="anode",size=ak.flatten(lctHMT[sel].WireNHits))
            #h.fill(sample="cathode",size=ak.flatten(lctHMT[sel].ComparatorNHits))
            
            h = hist.Hist("Events",hist.Bin("AnodeWireHit","AnodeWireHit",100,0,100),hist.Bin("CathodeHit","CathodeHit",100,0,100))
            sel = (lctHMT.sr==sr_map[i][0])|(lctHMT.sr==sr_map[i][1])
            h.fill(AnodeWireHit=ak.flatten(lctHMT[sel].WireNHits),
                   CathodeHit=ak.flatten(lctHMT[sel].ComparatorNHits))

            output[tag+"_"+label] = h

        return 

    def fillbasic(self,cls,dataset):
        version = "new"
        h1= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("nClusters", "nClusters", 10, 0, 10))
        h2= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("eta", "eta", 40, -5, 5))    
        h3= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("phi", "phi", 40, -np.pi, np.pi))        
        h4= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("size", "size", 100, 0, 1000))        
        h5= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("time", 'ClusterTime[ns]',80, -80, 80))            
        if version=="DT":
            h6 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("nMB1", 'nMB1',40, 0, 40))
            h12= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("nMB2", "nMB2", 40, 0, 40))        
        else:
            h6 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("nME11_12", 'nME11_12',40, 0, 40))        
            h12= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("timeSpread", "timeSpread[ns]", 80, 0, 100))
            
        h7 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("x", "x", 100, -1000, 1000))        
        h8 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("y", "y", 100, -1000, 1000))        
        h9 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("z", "z", 100, -1000, 1000))            
        h10= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("Nstation10", "Nstation10", 8, 0, 8))    
        h11= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("AvgStation10", "AvgStation10", 40, 0, 8))  

        hall = [h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12]
        h1.fill(sample=dataset,nClusters=ak.num(cls))
        h2.fill(sample=dataset,eta = ak.flatten(cls.Eta))
        h3.fill(sample=dataset,phi = ak.flatten(cls.Phi)) 
        h4.fill(sample=dataset,size = ak.flatten(cls.Size))
        h5.fill(sample=dataset,time = ak.flatten(cls.Time)) 
        if version=="DT":
            h6.fill (sample=dataset,nMB1 = ak.flatten(cls.nMB1))        
            h12.fill(sample=dataset,nMB2 = ak.flatten(cls.nMB2))                
        elif version=="new":
            h6.fill (sample=dataset,nME11_12 = ak.flatten(cls.ME11_12))
            h12.fill(sample=dataset,timeSpread = ak.flatten(cls.TimeSpread))                
    
        h7.fill( sample=dataset,x = ak.flatten(cls.X))    
        h8.fill( sample=dataset,y = ak.flatten(cls.Y))    
        h9.fill( sample=dataset,z = ak.flatten(cls.Z))        
        h10.fill(sample=dataset,Nstation10 = ak.flatten(cls.Nstation10))            
        h11.fill(sample=dataset,AvgStation10 = ak.flatten(cls.AvgStation10)) 

        return hall 


    def process(self, events):

        dataset = events.metadata['dataset']
        events = self.pack(events)
        cls = events.cls
        lctHMT = ak.zip({k.replace("lctHMT_",""):getattr(events,k) for k in events.fields if k.startswith("lctHMT_")})
        elctHMT = ak.zip({k.replace("elctHMT_",""):getattr(events,k) for k in events.fields if k.startswith("elctHMT_")})

        events.uniqueStation = ak.max(events.cscRechitsChamber,axis=1)==ak.min(events.cscRechitsChamber,axis=1)
        events.passL1_emul = ak.any(events.elctHMT_bits>1,axis=1)
        cls.R = (cls.X**2+ cls.Y**2)**0.5
        #mu_veto = ~(ak.all( (events.muonPt>10)&(events.muonIsGlobal)&(events.muonIsMedium) ,axis=1) )

        binning = np.array(
            [ 0.,  20.,  40.,  60.,  80., 100., 120., 140., 160., 180., 200.,
           220., 240., 260., 280., 300., 500., 700.]
        )

        ME11 = (cls.R>100)&(cls.R<275) &(abs(cls.Z)>580)&(abs(cls.Z)<632) 
        ME12 = (cls.R>275)&(cls.R<465) &(abs(cls.Z)>668)&(abs(cls.Z)<724)
        ME13 = (cls.R>505)&(cls.R<700) &(abs(cls.Z)>668)&(abs(cls.Z)<724)
        
        ME21 = (cls.R>139)&(cls.R<345) &(abs(cls.Z)>789)&(abs(cls.Z)<850)
        ME22 = (cls.R>357)&(cls.R<700) &(abs(cls.Z)>791)&(abs(cls.Z)<850)
        
        ME31 = (cls.R>160)&(cls.R<345) &(abs(cls.Z)>915)&(abs(cls.Z)<970)
        ME32 = (cls.R>357)&(cls.R<700) &(abs(cls.Z)>911)&(abs(cls.Z)<970)
        
        ME41 = (cls.R>178)&(cls.R<345) &(abs(cls.Z)>1002)&(abs(cls.Z)<1063)
        ME42 = (cls.R>357)&(cls.R<700) &(abs(cls.Z)>1002)&(abs(cls.Z)<1063)

        regions = {        "ME11":ME11,         "ME12":ME12,         "ME13":ME13, 
                            "ME21":ME21,         "ME22":ME22, 
                            "ME31":ME31,         "ME32":ME32, 
                            "ME41":ME41,         "ME42":ME42, 
        }

        all_thresholds = {
            "ME11":[120,100,12],
            "ME12":[120,100,12],
            "ME13":[7],
            "ME21":[23],
            "ME22":[12],
            "ME31":[21],
            "ME32":[12],
            "ME41":[25],
            "ME42":[12],

        }
        h= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("ClusterSize", "ClusterSize", binning)) 
        ## output hist
        #output = {}
        output = self.accumulator.identity()  ## get from histograms
        output["sumw"][dataset] += len(events)
        hall = self.fillbasic(cls,dataset)
        for h in hall:
            key = h.dense_axes()[0].name
            output[key] = h        

        self.fillHMT(elctHMT,output,"elctHMT")

        # default thresholds
        for regName,region in regions.items():
            sel = ak.any(region&(cls.Nstation10==1),axis=1) & (events.runNum>360019)
            denom = ak.fill_none(ak.max(cls[sel].Size,axis=1),-1)
            numer = ak.fill_none(denom[events[sel].passL1==1],-1)   
            output[regName] = self.addEff(numer,denom)
        # alternative thresholds
        for regName,thresholds in all_thresholds.items():
            region = regions[regName]
            sel = ak.all(region&(cls.Nstation10==1),axis=1) 
            denom = ak.fill_none(ak.max(cls[sel].Size,axis=1),-1)
            for threshold in thresholds:
                numer = ak.fill_none(denom[ak.any(elctHMT[sel].WireNHits>threshold,axis=1)],-1)   
                hname = "_".join([regName,str(threshold)])
                output[hname] = self.addEff(numer,denom)
 
        return output

    def addEff(self,numer,denom):
        h= hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("ClusterSize", "ClusterSize", 50, 10, 1000))        
        h.fill(sample="numer",ClusterSize=numer)
        h.fill(sample="denom",ClusterSize=denom)
        return h


    def postprocess(self, accumulator):
        pass
