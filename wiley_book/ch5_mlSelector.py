import simuOpt
simuOpt.setOptions(quiet=True, alleleType='binary')
import simuPOP as sim

pop = sim.Population(size=[10000]*3, loci=[1]*2, infoFields='fitness')
simu = sim.Simulator(pop, rep=3)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    preOps=[
        sim.MaSelector(loci=0, fitness=[1, 0.99, 0.98], reps=0),
        sim.MaSelector(loci=0, fitness=[1, 1, 0.99], reps=1),
        sim.MlSelector([
            sim.MaSelector(loci=0, fitness=[1, 0.99, 0.98]),
            sim.MaSelector(loci=1, fitness=[1, 1, 0.99])],
            mode=sim.MULTIPLICATIVE, reps=2)
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=[0,1], step=50),
        sim.PyEval(r"'%.3f\t' % alleleFreq[0][1]", reps=0, step=50),
        sim.PyEval(r"'%.3f\t' % alleleFreq[0][1]", reps=1, step=50),
        sim.PyEval(r"'%.3f\t%.3f\n' % (alleleFreq[0][1], alleleFreq[1][1])",
            reps=2, step=50),
    ],
    gen = 151
)

