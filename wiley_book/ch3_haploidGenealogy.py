import simuPOP as sim
pop = sim.Population(1000, ploidy=1, ancGen=-1,
            infoFields=['ind_id', 'father_id'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
    ],
    matingScheme=sim.RandomSelection(
        ops=[
            sim.IdTagger(),
            sim.PedigreeTagger(infoFields='father_id')
        ],
    ),
    gen = 1000
)
# a pedigree with only paternal information
pop.asPedigree(motherField='')
IDs = pop.identifyAncestors()
allIDs = [ind.ind_id for ind in pop.allIndividuals()]
removedIDs = list(set(allIDs) - set(IDs))
pop.removeIndividuals(IDs=removedIDs)
# number of ancestors...
sizes = [pop.popSize(ancGen=x) for x in range(pop.ancestralGens())]
print(sizes[0], sizes[100], sizes[500], sizes[999])
