import simuPOP as sim
def calcNe(pop, param):
    'Calculated effective number of disease alleles at specified loci (param)'
    sim.stat(pop, alleleFreq=param)
    ne = {}
    for loc in param:
        freq = pop.dvars().alleleFreq[loc]
        sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
        if sumFreq == 0:
            ne[loc] = 0
        else:
            ne[loc] = 1. / sum([(freq[x]/sumFreq)**2 \
                for x in list(freq.keys()) if x != 0])
    # save the result to the sim.Population.
    pop.dvars().ne = ne
    return True

pop = sim.Population(1000, loci=2)
sim.initGenotype(pop, freq=[0.9] + [0.01]*10, loci=0)
sim.initGenotype(pop, freq=[0.9] + [0.05] + [0.005]*10, loci=1)
calcNe(pop, param=[0,1])
print(pop.dvars().ne)
