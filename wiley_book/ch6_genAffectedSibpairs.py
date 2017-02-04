import simuPOP as sim

from simuPOP.sampling import drawAffectedSibpairSample

def genAffectedSibpairSample(pop, nFamilies, penetrance):
    '''Draw nFamilies affected sibpairs and their parents by producing
    siblings from pop repeatedly until enough affected sibpairs are
    collected. A penetrance operator is needed to assign affection status
    to each offspring.
    '''
    pop1 = pop.clone()
    pop1.setAncestralDepth(1)
    pop1.addInfoFields(['ind_id', 'father_id', 'mother_id'])
    pop1.evolve(
        initOps=sim.IdTagger(),
        matingScheme=sim.RandomMating(
            ops=[
                sim.MendelianGenoTransmitter(),
                penetrance,
                sim.IdTagger(),
                sim.PedigreeTagger(),
            ],
            numOffspring=2,
            subPopSize=pop.popSize()*2
        ),
        gen=1
    )
    sim.stat(pop1, numOfAffected=True)
    return drawAffectedSibpairSample(pop1, nFamilies)

if __name__ == '__main__':
    pop = sim.Population(size=10000, loci=1)
    sim.initGenotype(pop, freq=[0.5, 0.5])
    sim.initSex(pop)
    sim.maPenetrance(pop, loci=0, penetrance=[0.05, 0.15, 0.30])
    sample = genAffectedSibpairSample(pop, 100, 
        sim.MaPenetrance(loci=0, penetrance=[0.05, 0.15, 0.30]))
    #
    sim.stat(sample, numOfAffected=True)
    print(sample.dvars().numOfAffected, sample.dvars().numOfUnaffected)

