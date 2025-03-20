import importlib.resources 
import pickle

def build_lumimask(filename):
    from coffea.lumi_tools import LumiMask
    #print(filename)
    with importlib.resources.path("HMTprocessor.data", filename) as path:
        return LumiMask(path)
lumiMasks = {
    "2016": build_lumimask("Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
    "2017": build_lumimask("Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
    "2018": build_lumimask("Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
    "2022": build_lumimask("Cert_Collisions2022_355100_362760_Golden.txt"),
    "2023": build_lumimask("Cert_Collisions2023_366442_370790_Golden.txt"),
    #"2024": build_lumimask("Cert_Collisions2024_378981_385194_Golden.txt"),
    "2024": build_lumimask("Cert_Collisions2024_378981_386951_Golden.json"),
}

#def get_lumimask():
#    with importlib.resources.path("HMTprocessor.data", "lumi_masks.pkl") as path:
#        with open(path, "rb") as handle:
#            return pickle.load(handle)
#
#lumiMasks = get_lumimask()
