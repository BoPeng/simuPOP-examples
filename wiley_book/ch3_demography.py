import simuPOP as sim
from math import ceil
def demoModel(N0, N1, G0):
    def func(gen):
        if gen < G0:
            return N0
        else:
            return N1
    return func

def simulate(demo, gen):
    pop = sim.Population(size=demo(0))
    pop.evolve(
        initOps=sim.InitSex(),
        preOps=[
            sim.Stat(popSize=True),
            sim.PyEval(r"'%d: %d ' % (gen, popSize)"),
        ],
        matingScheme=sim.RandomMating(subPopSize=demo),
        postOps=[
            sim.Stat(popSize=True),
            sim.PyEval(r"'--> %d\n' % popSize"),
        ],
        gen=gen
    )

simulate(demoModel(100, 1000, 2), 5)

