from coffea import  processor
from optparse import OptionParser
import coffea
import pickle,glob
import time
from coffea.nanoevents import NanoEventsFactory, BaseSchema,NanoAODSchema
from HMTprocessor.HMTproc import MyProcessor
from HMTprocessor.MDSnano import MDSnanoProcessor
from HMTprocessor.MDSnanoSchema import MDSNanoAODSchema
#from HMTprocessor.MDStrigProcessor import MDStrigProcessor
#from HMTprocessor.BaseMDSnanoProcessor import BaseMDSnanoProcessor

def runLocal(outf="test.pickle",fileset="test.json",**options):
    #p = MDStrigProcessor()
    p = MyProcessor()

    if options['full']:      
        out = processor.run_uproot_job(
            fileset,
            treename='simpleCSCshowerTreeMaker/hmt',
            #treename='Events',
            processor_instance=p,
            executor=processor.iterative_executor,
            executor_args={
                "schema": BaseSchema,
                #"schema": MDSNanoAODSchema,
            },
            chunksize=10000,
        )
    else:
        print("loading 1 chunk")
        iterative_run = processor.Runner(
            executor = processor.IterativeExecutor(compression=None),
            schema=BaseSchema,
            #schema=MDSNanoAODSchema,
            maxchunks=1,
            chunksize=100,
        )
        out =iterative_run(
                fileset,
                treename='simpleCSCshowerTreeMaker/hmt',
                #treename='Events',
                processor_instance=p,
            )

    print(out)
    with open(outf,'wb') as f:
        pickle.dump(out,f)
    return

def runLPC(outf="test.pickle",fileset="test.json",**options):
    import time
    from distributed import Client
    from lpcjobqueue import LPCCondorCluster
    tic = time.time()
    #cluster = LPCCondorCluster(log_directory="/uscms/home/kkwok/log")
    #cluster = LPCCondorCluster()
    #cluster = LPCCondorCluster(shared_temp_directory="/tmp")
    #cluster = LPCCondorCluster(shared_temp_directory="/tmp", memory='6GB',
    cluster = LPCCondorCluster(shared_temp_directory='/tmp', 
                                 worker_extra_args=['--worker-port 10000:10070', '--nanny-port 10070:10100', '--no-dashboard'],
                                 job_script_prologue=[],
                                transfer_input_files=["HMTprocessor/"]) 

    # minimum > 0: https://github.com/CoffeaTeam/coffea/issues/465
    cluster.adapt(minimum=4, maximum=100)
    client = Client(cluster)

    exe_args = {
        "client": client,
        "savemetrics": True,
        "schema": BaseSchema,
        #"schema": MDSNanoAODSchema,
        "align_clusters": False,
    }

    #client.upload_file("HMTprocessor.zip")

    p = MyProcessor()
    #p = MDSnanoProcessor()
    #p = MDStrigProcessor()

    client.wait_for_workers(4)
    hists, metrics = processor.run_uproot_job(
        fileset,
        treename='simpleCSCshowerTreeMaker/hmt',
        #treename='Events',
        processor_instance=p,
        executor=processor.dask_executor,
        executor_args=exe_args,
        chunksize=10000,
    )


    elapsed = time.time() - tic
    print(f"Output: {hists}")
    print(f"Metrics: {metrics}")
    print(f"Finished in {elapsed:.1f}s")
    print(f"Events/s: {metrics['entries'] / elapsed:.0f}")

    with open(outf,'wb') as f:
        pickle.dump(hists,f)
    return

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('--test', dest='test', action='store_true',default = False, help='Run local test with small fileset')
    parser.add_option('--local', dest='local', action='store_true',default = False, help='Run local test with 1 chunk of full fileset')
    parser.add_option('--condor', dest='condor', action='store_true',default = False, help='Run local test with 1 chunk of full fileset')
    parser.add_option('--fileset', dest='fileset',default = "test.json", help='Run local test with 1 chunk of full fileset')
    parser.add_option('--full', dest='full', action='store_true',default = False, help='Run full file chunks')
    parser.add_option('-o', dest='outf', default='histograms.pickle', help='collection of histograms')

    (options, args) = parser.parse_args()
    procOptions       = vars(options)
    procOptions       = {k:v for k,v in procOptions.items() if k not in ["fileset","outf","isElectronChannel"]} ## write these 3 option explicitly

    outf    = options.outf
    print(" Coffea version = ", coffea.__version__)

    #fileset ="/eos/uscms/store/user/kkwok/llp/MDS/MuonRun3_MDSNtupler_V1p19_Data2022_Run2022E-PromptReco-v1_v1_v1_goodLumi.root"
    #fileset ="root://cmseos.fnal.gov//store/user/kkwok/llp/MDS/MuonRun3_MDSNtupler_V1p19_Data2022_Run2022E-PromptReco-v1_v1_v1_goodLumi.root"
    #fileset ="Run2022E.json"
    #fileset ="ZeroBias.json"
    #fileset ="Muon2024.json"
    fileset = options.fileset
 
    if options.test:
        runLocal("test.pickle","test.json",**procOptions)
    elif options.local:
        print("full               = ", options.full)
        runLocal(outf,fileset,**procOptions)
        pass
    elif options.condor:
        runLPC(outf,fileset,**procOptions)
        pass
