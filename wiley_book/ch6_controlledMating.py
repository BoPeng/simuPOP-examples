import simuPOP as sim
from simuPOP.utils import Trajectory, simulateBackwardTrajectory
from math import exp
def Nt(gen):
    'An exponential population expansion model'
    return int(5000 * exp(.00115 * gen))

# simulate a trajectory backward in time, from generation 1000
traj = simulateBackwardTrajectory(N=Nt, fitness=[1, 1.01, 1.02], nLoci=2,
     endGen=1000, endFreq=[0.1, 0.2])
# print out mutants in the format of (loc, gen, subPop)
print(traj.mutants())
pop = sim.Population(size=Nt(0), loci=[1]*2)
# save Trajectory function in the sim.population's local namespace
# so that the sim.PyEval operator can access it.
pop.dvars().traj = traj.func()
pop.evolve(
    initOps=[sim.InitSex()],
    preOps=traj.mutators(loci=[0, 1]),
    matingScheme=sim.ControlledRandomMating(loci=[0, 1], alleles=[1, 1],
        subPopSize=Nt, freqFunc=traj.func()),
    postOps=[
        sim.Stat(alleleFreq=[0, 1], begin=500, step=100),
        sim.PyEval(r"'%4d: %.3f (exp: %.3f), %.3f (exp: %.3f)\n' % (gen, alleleFreq[0][1],"
            "traj(gen)[0], alleleFreq[1][1], traj(gen)[1])",
            begin=500, step=100)
    ],
    gen=1001  # evolve 1001 generations to reach the end of generation 1000
)
