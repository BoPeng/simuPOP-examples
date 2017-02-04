import simuPOP as sim
pop = sim.Population(100, loci=[5]*3, infoFields='parent_idx')
pop.evolve(
    initOps=sim.InitGenotype(freq=[0.2]*5),
    preOps=sim.Dumper(structure=False, max=5),
    matingScheme=sim.HomoMating(
        sim.SequentialParentChooser(),
        sim.OffspringGenerator(ops=[
            sim.SelfingGenoTransmitter(),
            sim.ParentsTagger(infoFields='parent_idx'),
        ])
    ),
    postOps=sim.Dumper(structure=False, max=5),
    gen = 1
)
