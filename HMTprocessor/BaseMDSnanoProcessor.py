import awkward as ak
from coffea import processor
from coffea.nanoevents.methods import candidate
#import hist
from coffea import hist
import warnings
import numpy as np
from HMTprocessor import util
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods import vector
from HMTprocessor.corrections import lumiMasks
ak.behavior.update(candidate.behavior)

class BaseMDSnanoProcessor(processor.ProcessorABC):
    """
    Base Processor for MDS nano for  muon shower in the muon system. 

    """
    def prep(self,arr,fill_none=-1):
        return ak.flatten(ak.fill_none(arr,fill_none))

    def fillbasic(self,samples={},hall=None):
        if hall==None:
            h1= hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("nClusters", "nClusters", 10, 0, 10))
            h2= hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("eta", "eta", 40, -5, 5))    
            h3= hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("phi", "phi", 64, -3.2, 3.2))        
    
            bin_size = [   0.,  50.,   60.,   80.,  100.,  120.,  140.,  160.,
            180.,  200.,  220.,  240.,  260.,  280.,  300.,  320.,  340.,
            360.,  380.,  400.,  420.,  440.,  460.,  480.,  500.,  520.,
            540.,  560.,  580.,  600.,  620.,  640.,  660.,  680.,  700.,
            720.,  740.,  760.,  780.,  800.,  820.,  840.,  860.,  880.,
            900.,  920.,  940.,  960.,  980., 1000., 1100., 1200., 1300, 1400,1500,2000,4000.]
    
            h4=   hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("size", "size", bin_size))        
            h5=   hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("time", 'ClusterTime[ns]',80, -80, 80))            
            h6 =  hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("nME11_12", 'nME11_12',40, 0, 40))        
            h7 =  hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("x", "x", 50, -1000, 1000))        
            h8 =  hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("y", "y", 50, -1000, 1000))        
            h9 =  hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("z", "z", 50, -1000, 1000))            
            h10=  hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("Nstation10", "Nstation10", 8, 0, 8))    
            h11=  hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("AvgStation10", "AvgStation10", 40, 0, 8))  
            h12=  hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("timeSpread", "timeSpread[ns]", 50, 0, 100))
            h13 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Cat("reg","reg"),hist.Bin("R", "R", 50, 0, 1000))            
            h14 = hist.Hist("Events",hist.Cat("sample","sample"),hist.Bin("eta", "eta", 30,-3,3), hist.Bin("phi", "phi", 32,0,3.2))

            hall={"h1":h1,"h2":h2,"h3":h3,"h4":h4,"h5":h5,"h6":h6,"h7":h7,"h8":h8,"h9":h9,"h10":h10,"h11":h11,"h12":h12,"h13":h13,"h14":h14}

        for (dataset,regName),cls in samples.items():
            hall['h1'].fill(  sample=dataset,reg=regName,nClusters=ak.num(cls))
            hall['h2'].fill(  sample=dataset,reg=regName,eta = self.prep(cls.eta,-999))
            hall['h3'].fill(  sample=dataset,reg=regName,phi = self.prep(cls.phi,-999)) 
            hall['h4'].fill(  sample=dataset,reg=regName,size = self.prep(cls.size,))
            hall['h5'].fill(  sample=dataset,reg=regName,time = self.prep(cls.time,-999)) 
            hall['h6'].fill(  sample=dataset,reg=regName,nME11_12 = self.prep(cls.nME11+cls.nME12))
            hall['h7'].fill(  sample=dataset,reg=regName,x = self.prep(cls.x,-9999))    
            hall['h8'].fill(  sample=dataset,reg=regName,y = self.prep(cls.y,-9999))    
            hall['h9'].fill(  sample=dataset,reg=regName,z = self.prep(cls.z,-9999))
            hall['h10'].fill( sample=dataset,reg=regName,Nstation10 = self.prep(cls.nStation))            
            hall['h11'].fill( sample=dataset,reg=regName,AvgStation10 = self.prep(cls.avgStation)) 
            hall['h12'].fill( sample=dataset,reg=regName,timeSpread = self.prep(cls.timeSpread))                
            hall['h13'].fill( sample=dataset,reg=regName,R = self.prep((cls.x**2+cls.y**2)**0.5))
            hall['h14'].fill( sample=dataset,eta = self.prep(cls.eta,-999), phi = self.prep(cls.phi,-999))

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

    def compute_ctau(self,llp):
        decay_vx = llp.children.vx
        decay_vy = llp.children.vy
        decay_vz = llp.children.vz
    
        beta = np.sqrt(llp.p2)/llp.E
        decay_v = np.sqrt(decay_vx**2+decay_vy**2+decay_vz**2)
        gamma = 1.0/np.sqrt(1-beta**2)
        ctau = decay_v / (beta*gamma)    

        #reduce over childern's axis
        ctau = ak.firsts(ctau,axis=2)    
        
        return ctau 

    def build_cls(self,cls,llp=None):
        cls = ak.with_field(cls,(cls.x**2+cls.y**2)**0.5,"R")
        ## sort cls in event according to nHit
        idx   = ak.argsort(cls.size,ascending=False)
        cls   = cls[idx] 
        if llp is not None:
            cls = util.match(cls,llp,"matched_llp",0.4)
            cls = util.dr_min(cls,llp,"dr_llp")
        return cls

    def build_llp(self,events,llpId=9000006):
        if self.isSignal:
            gen = events.GenPart
            llp = gen[abs(gen.pdgId)==llpId]
            llp['ctau'] = compute_ctau(llp)
            idx   = ak.argsort(llp.ctau,ascending=False)
            llp   = llp[idx] 
            
            llp = ak.with_field(llp,ak.firsts(llp.children.vx,axis=2),"X")
            llp = ak.with_field(llp,ak.firsts(llp.children.vy,axis=2),"Y")
            llp = ak.with_field(llp,ak.firsts(llp.children.vz,axis=2),"Z")
            llp = ak.with_field(llp,np.sqrt(llp.X**2+llp.Y**2),"R")

            csc = (abs(llp.eta)<2.4) &(abs(llp.eta)>0.9)& (llp.R<695.5) &( abs(llp.Z)>400) &( abs(llp.Z)<1100)
            dt  = (abs(llp.Z)<661) & (llp.R>200) & (llp.R<800)
            llp['csc'] = csc
            llp['dt'] = dt  

            return llp
        return None 
       
    def AddNoiseFilters(self,events,Usefilters=[]):
        flags = events.Flag
        #From https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Run_3_recommendations
        allfilters = [
            "goodVertices_pRECO",
            "globalSuperTightHalo2016Filter_pRECO",
            "EcalDeadCellTriggerPrimitiveFilter_pRECO",
            "BadPFMuonDzFilter_pRECO",
            "BadPFMuonFilter_pRECO",
            "hfNoisyHitsFilter_pRECO",
            "eeBadScFilter_pRECO",
        ]
        if len(Usefilters)>0:
            allUseFilters = np.array([ getattr(flags,f) for f in allfilters if f in Usefilters])
        else:
            allUseFilters = np.array([ getattr(flags,f) for f in allfilters ])

        # transposed axis = flag dim, ask all flags to be true
        return ak.all(allUseFilters.transpose(),axis=1)

    def ApplyGoldenJSON(self,events): 
        if self.isData:
            dataset = events.metadata['dataset']
            year = None
            if "_2016" in dataset: year = "2016"
            elif "_2017" in dataset: year = "2017"
            elif "_2018" in dataset: year = "2018"
            elif "2022" in dataset: year = "2022"
            elif "2023" in dataset: year = "2023"
            elif "2024" in dataset: year = "2024"
            else:
                warnings.warn(" %s is data, but does not contain one of the strings: [_2016,_2017,_2018, 2022, 2023, 2024]. No golden json mask applied." % dataset, RuntimeWarning)
            self._year = year
            if year is not None:
                events = events[lumiMasks[year](events.run,events.luminosityBlock)]
            return events
        else:
            pass

    def __init__(self):
        self.isSignal = False
        self.isData   = True
        
    @property
    def accumulator(self):
        return self._accumulator
