import simuPOP as sim
def genCaseControlSample(pop, nCase, nControl, penetrance):
    '''Draw nCase affected and nControl unaffected individuals by producing
    offspring from pop repeatedly until enough cases and controls are
    collected. A penetrance operator is needed to assign affection status
    to each offspring.
    '''
    sample = pop.clone()
    sample.setVirtualSplitter(sim.ProductSplitter([
        sim.AffectionSplitter(),
        sim.RangeSplitter([[0, nCase], [nCase, nCase + nControl]])]))
    sample.evolve(
        matingScheme=sim.RandomMating(ops=[
            sim.MendelianGenoTransmitter(),
            penetrance,
            sim.DiscardIf(True, subPops=[(0,0), (0,3)])],
            subPopSize=nCase + nControl
        ),
        gen=1
    )
    return sample

if __name__ == '__main__':
    pop = sim.Population(size=10000, loci=1)
    sim.initGenotype(pop, freq=[0.8, 0.2])
    sim.initSex(pop)
    sample = genCaseControlSample(pop, 500, 500, 
        sim.MaPenetrance(loci=0, penetrance=[0.01, 0.02, 0.10]))
    #
    sim.stat(sample, numOfAffected=True)
    print(sample.dvars().numOfAffected, sample.dvars().numOfUnaffected)

