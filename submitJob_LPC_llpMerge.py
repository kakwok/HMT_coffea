import glob
import sys,  os 
from optparse import OptionParser,OptionGroup

def exec_me(command, dryRun=False):
    print(command)
    if not dryRun:
        os.system(command)

def write_condor(njobs , exe='runjob', files = [], dryRun=True):
    fname = '%s.jdl' % exe
    out = """universe = vanilla
Executable = {exe}.sh
Should_Transfer_Files = YES 
WhenToTransferOutput = ON_EXIT_OR_EVICT
request_memory = 4GB
Transfer_Input_Files = {exe}.sh,{files}
Output = {exe}.$(Process).$(Cluster).stdout
Error  = {exe}.$(Process).$(Cluster).stdout
Log    = {exe}.$(Process).$(Cluster).log
Arguments = $(Process) {njobs}
Queue {njobs}
    """.format(exe=exe, files=','.join(files), njobs=njobs)
    with open(fname, 'w') as f:
        f.write(out)
    if not dryRun:
        os.system("condor_submit %s" % fname)


def write_bash(temp = 'runjob.sh', eoscp="" ):
    out = '#!/bin/bash\n'
    out += 'date\n'
    out += 'MAINDIR=`pwd`\n'
    out += 'echo "Process = $1"\n'
    out += 'echo "njobs = $2"\n'
    out += 'source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc12-opt/setup.sh\n'
    out += 'ls\n'
    out += 'python ${MAINDIR}/hadd_group.py -l ${MAINDIR}/input_list.txt -n $2 -i $1 -o merged_$1.root\n'
    if eoscp!="":
        out += 'echo "coping to eos: "+%s  \n'%eoscp
        out +=  eoscp + '\n'
    out += 'echo "Inside $MAINDIR:"\n'
    out += 'ls\n'
    if eoscp!="":
        out += 'cd $MAINDIR  \n'
        out += 'echo "remove output local file"  \n'
        out += 'rm -rf *.root \n'
    out += 'ls\n'
    out += 'date\n'
    with open(temp, 'w') as f:
        f.write(out)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--clean', dest='clean', action='store_true',default = False, help='clean submission files', metavar='clean')
    parser.add_option('--dryRun', dest='dryRun', action='store_true',default = False, help='write submission files only', metavar='dryRUn')
    parser.add_option('-o', '--odir', dest='odir', default='./', help='directory to write histograms/job output', metavar='odir')
    parser.add_option('-d', '--dataset', dest='dataset', default='EGamma0')
    parser.add_option('--era',  dest = 'era',default='MDSnano_Run2024B-PromptReco-v1')

    (options, args) = parser.parse_args()
    dryRun= options.dryRun 

    dataset  = options.dataset
    era      = options.era
    #eras     = glob.glob(f"/uscms/home/kkwok/lpclonglived/MDSnano/{dataset}/*")
    #eraPath  = f"/uscms/home/kkwok/lpclonglived/MDSnano/{dataset}/{era}" 
    eras     = glob.glob(f"/uscms/home/kkwok/lpclonglived/HLT/zerobias24/{dataset}/*")
    eraPath  = f"/uscms/home/kkwok/lpclonglived/HLT/zerobias24/{dataset}/{era}" 
    disk={}

    # prepare in main dir
    with open("./disk.txt","r") as file:
        for line in file:
            size,key = line.split()[0],line.split()[1]
            disk[key] = size
    #[ print(era, disk[era]) for era in eras]

    nGB            = int(disk[eraPath].replace("G",""))
    with open("input_list.txt",'w') as fout:
        #files = glob.glob(f"/uscms/home/kkwok/lpclonglived/MDSnano/{dataset}/{era}/*/*/reco*.root")
        files = glob.glob(f"/uscms/home/kkwok/lpclonglived/HLT/zerobias24/{dataset}/{era}/*/*/HMT*.root")
        nfiles = len(files)
        for line in files: 
            eosline = line.replace("/uscms/home/kkwok/","root://cmseos.fnal.gov//store/user/")
            fout.write(eosline+"\n")

    nFilesPerChunk = int(nfiles/nGB)  # target ~ 1GB per file
    print(f"nfile = {nfiles}, Size = {nGB} GB, nFilesPerChunk = {nFilesPerChunk}")

    ##Small files used by the exe
    transfer_files = [
        os.getcwd()+"/hadd_group.py",
        "./input_list.txt",
    ]

    #submit from sub-dir
    exec_me(f"mkdir -p ./clean_up/{dataset}/{era}", False)
    exec_me(f"mv input_list.txt ./clean_up/{dataset}/{era}",False)
    os.chdir(f"./clean_up/{dataset}/{era}")

    eoscp        = 'xrdcp -f ${MAINDIR}/merged_$1.root %s'%(eraPath.replace("/uscms/home/kkwok/","root://cmseos.fnal.gov//store/user/"))

    exe = "runjob"
    write_bash(exe+".sh",eoscp=eoscp )
    write_condor(nGB, exe,  transfer_files, dryRun)
