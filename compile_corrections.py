
import pickle
from coffea.lumi_tools import LumiMask
if __name__ == "__main__":
    lumi_masks = {
        "2016": LumiMask("./HMTprocessor/data/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        "2017": LumiMask("./HMTprocessor/data/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
        "2018": LumiMask("./HMTprocessor/data/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
        "2022": LumiMask("./HMTprocessor/data/Cert_Collisions2022_355100_362760_Golden.txt"),
        "2023": LumiMask("./HMTprocessor/data/Cert_Collisions2023_366442_370790_Golden.txt"),
        #"2024": LumiMask("./HMTprocessor/data/Cert_Collisions2024_378981_385194_Golden.txt"),
        "2024": LumiMask("./HMTprocessor/data/Cert_Collisions2024_378981_386951_Golden.json"),
    }
    with open("lumi_masks.pkl", "wb") as handle:
        pickle.dump(lumi_masks, handle, protocol=pickle.HIGHEST_PROTOCOL)
