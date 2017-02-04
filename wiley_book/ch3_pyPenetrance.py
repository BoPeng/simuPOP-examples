import simuPOP as sim
import random

def breatCancer(geno, age):
    # brca1 is recessive
    if 0 in geno:
        return 0
    return [0, 0, 0.036, 0.18, 0.57, 0.75, 0.83][int(age/10)]

pop = sim.Population(size=10000, loci=1, infoFields='age',
    lociNames='BRCA1')
sim.initGenotype(pop, freq=[0.73, 0.27])
sim.initInfo(pop, lambda: random.randint(0, 69), infoFields='age')
sim.pyPenetrance(pop, func=breatCancer, loci=0)
pop.setVirtualSplitter(sim.CombinedSplitter([
    sim.InfoSplitter(field='age', cutoff=[30,40,50,60]),
    sim.ProductSplitter([
        sim.InfoSplitter(field='age', cutoff=[30,40,50,60]),
        sim.GenotypeSplitter(loci=0, alleles=[1,1], names='carrier'),
    ])]
))
sim.stat(pop, numOfAffected=True, subPops=[(0, sim.ALL_AVAIL)],
    vars=['propOfAffected', 'propOfAffected_sp'])
print('Population prevalence is %.2f%%' % (pop.dvars().propOfAffected*100))
for x in range(pop.numVirtualSubPop()):
    print('Prevalence in group %s is %.2f%%' % \
        (pop.subPopName((0,x)), pop.dvars((0,x)).propOfAffected*100))

