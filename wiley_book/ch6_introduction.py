import simuOpt
simuOpt.setOptions(quiet=True, alleleType='long')
import simuPOP as sim
pop = sim.Population(size=1000, loci=[1])
simu = sim.Simulator(pop, 10)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.PointMutator(loci=0, allele=1, inds=0)
    ],
    matingScheme=sim.RandomMating(),
    finalOps=sim.Stat(alleleFreq=0),
    gen = 100
)
print([x.dvars().alleleNum[0][1] for x in simu.populations()])
