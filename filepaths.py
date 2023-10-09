import glob,os
import json 

fileset = {
       'Run2022E':glob.glob("/eos/uscms/store/user/jongho/llp/pp/ZeroBias/ppZeroBias_cathode_hltzb_path_Run2022E-v1/230506_045208/000*/plots_*")        
}


def writejson(fileset,fout):
    outf = open(os.path.expandvars(fout),"w")
    finaljson={}
    for sample, filelist in fileset.items():
        print("working on ",sample)
        newpaths = []
        for f in filelist:
            #newpaths.append(f.replace("/eos/uscms/",'root://xcache//')) ## coffea-casa
            newpaths.append(f.replace("/eos/uscms/",'root://cmseos.fnal.gov//'))
        finaljson[sample] = newpaths
    outf.write((json.dumps(finaljson,indent=4)))
    print("Json written to :", fout)
    return

# For coffea-casa
def replaceFilepath(inputjson):
    outputjson = inputjson.replace(".json","_casa.json")
    with open(inputjson) as f:
        d = json.load(f)
    output ={}
    for k,flist in d.items():
        output[k]=[f.replace("cmsxrootd.fnal.gov","xcache") for f in flist]
    with open(outputjson, 'w') as outfile:
        outfile.write(json.dumps(output,indent=4))

writejson(fileset,"Run2022E.json")
