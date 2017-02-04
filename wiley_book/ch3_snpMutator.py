import simuPOP as sim
from math import log
pop = sim.Population(size=10000, loci=1)
simu = sim.Simulator(pop, rep=3)
simu.evolve(
    initOps=[
        sim.InitSex(),
        # Because of the use of parameter reps, these operators are
        # applied to different populations.
        sim.InitGenotype(freq=(0.2, 0.8), reps=0),
        sim.InitGenotype(freq=(0.5, 0.5), reps=1),
        sim.InitGenotype(freq=(0.8, 0.2), reps=2)],
    preOps=sim.SNPMutator(u=0.01, v=0.001),
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=0, step=100),
        sim.PyEval('gen', reps=0, step=100),
        sim.PyEval(r"'\t%.3f' % alleleFreq[0][0]", step=100),
        sim.PyOutput('\n', reps=-1, step=100)
    ],
    gen=500
)
