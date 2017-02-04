import simuPOP as sim
from random import randint

def podParentsChooser(pop, subPop):
    '''Choose parents of parents from different pods'''
    males = [x for x in pop.individuals(subPop) if x.sex() == sim.MALE]
    females = [x for x in pop.individuals(subPop) if x.sex() == sim.FEMALE]
    while True:
        # randomly choose a male
        male = males[random.randint(0, len(males)-1)]
        pod = male.pod
        # randomly choose a female from different pod
        while True:
            female = females[randint(0, len(females)-1)]
            if female.pod != pod:
                break
        yield (male, female)

pop = sim.Population(5000, loci=[1,1], infoFields=['pod'],
    chromTypes=[sim.AUTOSOME, sim.CUSTOMIZED])
pop.setVirtualSplitter(sim.InfoSplitter('pod', values=range(5)))
pop.evolve(
    initOps = [
        sim.InitSex(),
        # assign individuals to a random pod
        sim.InitInfo(lambda : randint(0, 4), infoFields='pod'),
        # only the first pod has the disease alleles
        sim.InitGenotype(freq=[0.8, 0.2], subPops=[(0,0)]),
    ],
    matingScheme = sim.HomoMating(
        sim.PyParentsChooser(podParentsChooser),
        sim.OffspringGenerator(numOffspring=1, ops=[
            sim.MendelianGenoTransmitter(),
            sim.MitochondrialGenoTransmitter(),
            # offspring stays with their natal pod
            sim.InheritTagger(mode=sim.MATERNAL, infoFields='pod')])),
    postOps = [
        # calulate allele frequency at each pod
        sim.Stat(alleleFreq=(0,1), vars='alleleFreq_sp',
            subPops=[(0, sim.ALL_AVAIL)]),
        sim.PyEval(r"'Loc0: %s Loc1: %s\n' % ("
    "', '.join(['%.3f' % subPop[(0,x)]['alleleFreq'][0][1] for x in range(5)]),"
    "', '.join(['%.3f' % subPop[(0,x)]['alleleFreq'][1][1] for x in range(5)]))"),
    ],
    gen = 10
)
