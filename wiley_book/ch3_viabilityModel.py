import simuPOP as sim
pop = sim.Population(size=10000, loci=2, infoFields='fitness')
a, b, c, d = 1, 1.5, 2.5, 4
r = 0.02
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=(0.5, 0.5)),
        sim.PyOutput('LD,   AB,  Ab,  aB,  ab,  avg fitness\n'),
    ],
    preOps=[
        sim.MaSelector(loci=[0, 1], wildtype=0,
            fitness=[a, b, a, c, d, c, a, b, a]),
        sim.Stat(meanOfInfo='fitness'),
    ],
    matingScheme=sim.RandomMating(ops=sim.Recombinator(rates=r)),
    postOps=[
        sim.Stat(haploFreq=[0,1], LD=[0,1], step=20),
        sim.PyEval(r"'%.3f %.2f %.2f %.2f %.2f %.2f\n' % (LD[0][1], "
            "haploFreq[(0,1)][(0,0)], haploFreq[(0,1)][(0,1)],"
            "haploFreq[(0,1)][(1,0)], haploFreq[(0,1)][(1,1)],"
            "meanOfInfo['fitness'])", step=20),
    ],
    gen=100
)
