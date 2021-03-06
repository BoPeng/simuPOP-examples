import simuPOP as sim
from random import randint
pop = sim.Population(10000, loci=1, 
    infoFields=['age', 'ind_id', 'father_id', 'mother_id'])
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[20, 50, 75]))
pop.evolve(
    initOps=[
        sim.InitSex(),
        # random assign age
        sim.InitInfo(lambda: randint(0, 74), infoFields='age'),
        sim.InitGenotype(freq=[0.5, 0.5]),
        # assign an unique ID to everyone.
        sim.IdTagger(),
    ],
    # increase the age of everyone by 1 before mating.
    preOps=sim.InfoExec('age += 1'),
    matingScheme=sim.HeteroMating([
        # all individuals with age < 75 will be kept. Note that
        # CloneMating will keep individual sex, affection status and all
        # information fields (by default).
        sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
        # only individuals with age between 20 and 50 will mate and produce
        # offspring. The age of offspring will be zero.
        sim.RandomMating(ops=[
            sim.IdTagger(),                   # give new born an ID
            sim.PedigreeTagger(),             # track parents of each individual
            sim.MendelianGenoTransmitter()],  # transmit genotype
            numOffspring=(sim.UNIFORM_DISTRIBUTION, 1, 3),
            subPops=[(0,1)])
    ]),
    gen = 200
)

from simuPOP import sampling
sample = sampling.drawNuclearFamilySample(pop, families=1, numOffspring=(2,3))
sim.dump(sample, structure=False)
