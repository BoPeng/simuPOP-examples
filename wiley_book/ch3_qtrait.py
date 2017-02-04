import simuPOP as sim
import random
pop = sim.Population(size=5000, loci=2, infoFields=['qtrait1', 'qtrait2', 'age'])
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[40]))
def qtrait(geno, age):
    'Return two traits that depends on genotype and age'
    return random.normalvariate(age * sum(geno), 10), random.randint(0, 10*sum(geno))

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.2, 0.8]),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        # use random age for simplicity
        sim.InitInfo(lambda:random.randint(20, 75), infoFields='age'),
        sim.PyQuanTrait(loci=(0,1), func=qtrait, infoFields=['qtrait1', 'qtrait2']),
        sim.Stat(meanOfInfo=['qtrait1'], subPops=[(0, sim.ALL_AVAIL)],
            vars='meanOfInfo_sp'),
        sim.PyEval(r"'Mean of trait1: %.3f (age < 40), %.3f (age >=40)\n' % "
            "(subPop[(0,0)]['meanOfInfo']['qtrait1'], subPop[(0,1)]['meanOfInfo']['qtrait1'])"),
    ],
    gen = 5
)

