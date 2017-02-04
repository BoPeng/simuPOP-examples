import simuPOP as sim
pop = sim.Population(size=10000, loci=10, lociPos=range(5) + range(10, 15))
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(haplotypes=[[0]*10,[1]*10]),
    ],
    matingScheme=sim.RandomMating(ops=sim.Recombinator(intensity=0.0005)),
    postOps=[
        sim.Stat(LD=[[1,2],[4,5],[8,9],[0,9]], step=10),
        sim.PyEval(r"'gen=%d\tLD12=%.3f (%.3f)\tLD45=%.3f (%.3f)\tLD09=%.3f\n'%"
            "(gen, LD[1][2], 0.25*0.9995**(gen+1), LD[4][5],"
            "0.25*0.9975**(gen+1),LD[0][9])", step=10)
    ],
    gen=100
)
