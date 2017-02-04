import simuPOP as sim

from simuPOP.sampling import drawCaseControlSample
# create a population with affected individuals
pop = sim.Population(size=10000, loci=5)
sim.initGenotype(pop, freq=[0.7, 0.3])
sim.maPenetrance(pop, loci=2, penetrance=[0.1, 0.12, 0.20])
sim.stat(pop, numOfAffected=True, association=range(5))
print(pop.dvars().numOfAffected)
# draw a case control sample and run association test
sample = drawCaseControlSample(pop, cases=500, controls=500)
sim.stat(sample, numOfAffected=True, association=range(5))
print(sample.dvars().numOfAffected, sample.dvars().numOfUnaffected)
print(', '.join(['%.3e' % sample.dvars().Allele_ChiSq_p[x] for x in range(5)]))
