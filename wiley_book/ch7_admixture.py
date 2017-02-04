import simuOpt
simuOpt.setOptions(gui=False, alleleType='binary')
import simuPOP as sim
pop.addInfoFields(['ancestry', 'migrate_to'])
# initialize ancestry
sim.initInfo(pop, [0]*pop.subPopSize(0) + [1]*pop.subPopSize(1),
    infoFields='ancestry')
# define two virtual subpopulations by ancestry value
pop.setVirtualSplitter(sim.InfoSplitter(field='ancestry', cutoff = [0.5]))
transmitters=[
    sim.MendelianGenoTransmitter(),
    sim.InheritTagger(mode=sim.MEAN, infoFields='ancestry')]
pop.evolve(
    initOps=sim.InitSex(),
    preOps=sim.Migrator(rate=[
        [0., 0], [0.05, 0]]), 
    matingScheme=sim.HeteroMating(
        matingSchemes=[
            sim.RandomMating(ops=transmitters),
            sim.RandomMating(subPops=[(0,0)], weight=-0.80, ops=transmitters),
            sim.RandomMating(subPops=[(0,1)], weight=-0.80, ops=transmitters)
        ],
    ),
    gen=10,
)
# remove the second subpop
pop.removeSubPops(1)
