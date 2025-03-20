#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import subprocess

import uproot3 as uproot
#import uproot 
import numpy as np
from coffea import lumi_tools

def getSamples(sampleInput):
    samples = {}
    #for era in ["A","B","C","D"]:
    for era,cmd in sampleInput.items():
        lines  = subprocess.getoutput(
                cmd
              #'ls /uscms/home/kkwok/lpclonglived/HNL/skim/SingleMuon_2018D/HeavyNeutralLepton_Tree_*.root'
              #'find s /eos/uscms/store/user/lpcbacon/15/JetHTRun2018%s*/*.root -not -empty -ls '%era
              ).split('\n')
        flist = [l.split()[-1] for l in lines]
        if 'directory' in flist:
            flist.remove('directory')
        samples[era]=flist
    return samples

def getSamplesFromList(inputjson):
    samples={}
    print("Calculating lumi for inputfile:", inputjson)
    with open(inputjson) as f:
        samples = json.load(f)
    for key,flist in samples.items():
        print(f'era = {key} , Nfile = ',len(flist))
    return samples

def get_lumilist(dataset, filename,treename="MuonSystem"):
     #print(filename)
     Badfiles = []
     #try:
     #   file = uproot.open(filename)
     #   if treename not in file.keys():
     #        print("Bad file:", filename)
     #        Badfiles.append(filename)
     #        return dataset, lumi_tools.LumiList(),Badfiles
     #except:
     #   print("Bad file:",filename)
     #   Badfiles.append(filename)
     #   return dataset, lumi_tools.LumiList(),Badfiles
     file = uproot.open(filename)

     if treename=="Events":
         #run, lumi = file["LuminosityBlocks"].arrays()['run'], file["LuminosityBlocks"].arrays()["luminosityBlock"]
         run, lumi = file["LuminosityBlocks"].array("run"), file["LuminosityBlocks"].array("luminosityBlock") 
     elif treename =="MuonSystem":
         tree = file[treename]
         run, lumi = tree["runNum"].array(), tree["lumiSec"].array()
     else:
         tree = file[treename]
         run, lumi = tree["runNum"].array(), tree["lumiNum"].array()
     if len(run)==0:
         Badfiles.append(filename)
         print("empty file:",filename)
         return dataset, lumi_tools.LumiList(),Badfiles
     lumilist = lumi_tools.LumiList(run, lumi)
     return dataset, lumilist, Badfiles

def printLumi(samples,treename="MuonSystem"):
    if np.all([("2016" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("metadata/lumi2016_new.csv")
    elif np.all([("2017" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("metadata/lumi2017.csv.gz")
    elif np.all([("2018" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("metadata/lumi2018.csv")
    elif np.all([("2022" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("metadata/lumi2022.csv")
    elif np.all([("2023" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("metadata/lumi2023.csv")
    elif np.all([("2024" in k) for k in samples.keys()]):
        lumivalues = lumi_tools.LumiData("./HMTprocessor/data/lumi2024_golden.csv")
    else:
        print(samples.keys()," must contain 2016,2017,2018,2022,2023 or 2024")
        return
    runs  = lumivalues._lumidata[:, 0].astype('u4')
    lumis = lumivalues._lumidata[:, 1].astype('u4')
    ll = lumi_tools.LumiList(runs,lumis)
    print(lumivalues.get_lumi(ll))

    tic = time.time()
    dataset_lumi = {}
    dataset_badFiles = {}
    nworkers = 8
    fileslice = slice(None)
    with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
        futures = set()
        for dataset, files in samples.items():
            futures.update(executor.submit(get_lumilist, dataset, file, treename) for file in files[fileslice])
        try:
            total = len(futures)
            processed = 0
            while len(futures) > 0:
                finished = set(job for job in futures if job.done())
                for job in finished:
                    dataset, accumulator ,badfiles = job.result()
                    if dataset in dataset_lumi:
                        dataset_lumi[dataset] += accumulator
                        dataset_badFiles[dataset] += badfiles
                    else:
                        dataset_lumi[dataset] = accumulator
                        dataset_badFiles[dataset] = badfiles
                    processed += 1
                    if processed % 100 == 0:
                        print("Processing: done with % 4d / % 4d files" % (processed, total))
                futures -= finished
            del finished
        except KeyboardInterrupt:
            print("Ok quitter")
            for job in futures: job.cancel()
        except:
            for job in futures: job.cancel()
            raise
    
    elapsed = time.time() - tic
    print(f"Finished in {elapsed:.1f}s")
   
    fname =[k for k in samples.keys() ][0]
    fout = open("lumi_%s.txt"%fname,"w") 
    print("dataset, lumi [/pb], lumisections, unique lumisections")
    fout.write("dataset, lumi [/pb], lumisections, unique lumisections\n")
    
    s = 0
    for ds, ll in dataset_lumi.items():
        lumi = lumivalues.get_lumi(ll.array)
        s+=lumi
        nunique = np.unique(ll.array, axis=0).shape[0]
        ntot = ll.array.shape[0]
        print("%s %0.2f %6d %6d" % (ds, lumi, ntot, nunique))
        fout.write("%s %0.2f %6d %6d\n" % (ds, lumi, ntot, nunique))
    print("total lumi ",s)
    fout.write("total lumi %s\n"%s)
    print("Bad files:  ")
    for ds, badfiles in dataset_badFiles.items():
        for f in badfiles:
            print(f)
            fout.write("%s\n" % (f))
 
   
if __name__ == '__main__':

    Treename = "Events"
    #printLumi(getSamplesFromList("test.json"),Treename)
    #printLumi(getSamplesFromList("Muon2024_Mar12.json"),Treename)
    printLumi(getSamplesFromList("Muon2024_Feb18.json"),Treename)
