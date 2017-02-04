import simuOpt
simuOpt.setOptions(gui=False, alleleType='binary')
import simuPOP as sim

import math
def penetrance(alpha, beta):
    def func(geno):
        g = geno[0] + geno[1]
        logit = alpha + beta*g
        return math.exp(logit) / (1. + math.exp(logit))
    return func

class TrioSampler:
    def __init__(self):
        # IDs of the parents of selected offspring
        self.parentalIDs = set()
    
    def _discardTrio(self, off):
        'Determine if the offspring can be kept.'
        if off.affected() and off.father_id not in self.parentalIDs and \
            off.mother_id not in self.parentalIDs:
            self.parentalIDs |= set([off.father_id, off.mother_id])
            return False
        # discard unaffected individual or individuals with duplicate parents
        return True
    
    def drawSample(self, pop, penet, nFamilies):
        self.pop = pop.clone()
        self.pop.addInfoFields(['ind_id', 'father_id', 'mother_id'])
        self.pop.setAncestralDepth(1)
        sim.tagID(self.pop, reset=True)
        self.pop.evolve(
            preOps = penet,
            matingScheme=sim.RandomMating(ops=[
                sim.MendelianGenoTransmitter(), # pass genotype
                sim.IdTagger(),       # assign new ID to offspring
                sim.PedigreeTagger(), # record the parent of each offspring
                penet,                # determine offspring affection status
                sim.DiscardIf(cond=self._discardTrio)
                ], subPopSize=nFamilies),
            gen = 1
        )
        return self.pop

sample = TrioSampler().drawSample(pop, 
    sim.PyPenetrance(penetrance(-0.5, -1.0), loci='rs2173746'), 1000)
