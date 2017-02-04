import simuOpt
simuOpt.setOptions(quiet=True, alleleType='long')
import simuPOP as sim
pop = sim.Population(size=1000, loci=[1])
simu = sim.Simulator(pop, 5)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.PyExec('introGen=[]')
    ],
    preOps=[
        sim.Stat(alleleFreq=0),
        sim.IfElse('alleleFreq[0][1] == 0', ifOps=[
            sim.PointMutator(loci=0, allele=1, inds=0),
            sim.PyExec('introGen.append(gen)')
        ]),
        sim.TerminateIf('alleleFreq[0][1] >= 0.05')
    ],
    matingScheme=sim.RandomMating()
)
# number of attempts
print([len(x.dvars().introGen) for x in simu.populations()])
# age of mutant
print([x.dvars().gen - x.dvars().introGen[-1] for x in simu.populations()])
