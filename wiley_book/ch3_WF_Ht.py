import simuPOP as sim
# [100]*100 means a population with 100 subpopulations, each of size 100.
pop = sim.Population([100]*100, loci=1)
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5]),
        sim.PyOutput('gen:   mean freq   mean Ht (expected Ht)\n')
        ],
    preOps=[
        # Statistics in subpopulations are by default not calculated.
        # This is changed by specifying variables with a '_sp' suffix.
        sim.Stat(alleleFreq=0, heteroFreq=0,
            vars=['alleleFreq_sp', 'heteroFreq_sp'], step=20),
        # Variables in subpopulations are stored in dictionaries
        # subPop[subPopID] where subPopID can be virtual.
        sim.PyEval(r'"%2d:    %.4f      %.4f (%.4f)\n" % (gen,'
            'sum([subPop[x]["alleleFreq"][0][1] for x in range(100)])/100.,'
            'sum([subPop[x]["heteroFreq"][0] for x in range(100)])/100.,'
            '0.5*(1-1/200.)**gen)', step=20)
    ],
    matingScheme=sim.RandomMating(),
    gen=100
)
