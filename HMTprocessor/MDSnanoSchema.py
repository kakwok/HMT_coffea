import awkward as ak
from coffea.nanoevents.methods import base, candidate, vector
from coffea.nanoevents.methods import nanoaod
from coffea.nanoevents import  NanoAODSchema
behavior = {}
behavior.update(nanoaod.behavior)

@ak.mixin_class(behavior)
class MDSClusterCollection(vector.LorentzVector, base.NanoCollection):
    @property
    def mass(self):
        return 0.0 * self.x
    
    @property
    def energy(self):
        return 0.0 * self.x

@ak.mixin_class(behavior)
class cscRechitsCollection(vector.LorentzVector, base.NanoCollection):
    @property
    def mass(self):
        return 0.0 * self.x
    
    @property
    def energy(self):
        return 0.0 * self.x
    

class MDSNanoAODSchema(NanoAODSchema):
    """MDSNano schema builder

    """
    mixins = {
        **NanoAODSchema.mixins,
        "MDSHLTCluster": "MDSClusterCollection",        
        "cscRechits": "cscRechitsCollection",        
        "cscMDSHLTCluster": "MDSClusterCollection",        
        "dtMDSHLTCluster": "MDSClusterCollection",        
    }
    all_cross_references = {
        **NanoAODSchema.all_cross_references,
    }
    @property
    def behavior(self):
        return behavior
