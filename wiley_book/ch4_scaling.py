def simuAscerBias(t, t0, N, N1, N2, v0, v1, v2, incProb, thresh):
    '''This function evolves a founder population of N individuals for
    t-t0 generations. A copy of this population is expanded to a size of
    N1 and continued to evolve for t-t0 generations. If a randomly chosen
    allele from the this population has more tandem repeats than a
    specified threshold, the same copy of population (at generation
    t-t0) is evolved similarly with a population size of N1. Ascertainment
    bias is returned as the length difference between the first and 
    second populations.'''
    while True:
        # Evolve an ancestral population with size N0
        pop = sim.Population(size=N, loci=1)
        pop.evolve(
            initOps=[
                sim.InitSex(),
                sim.InitGenotype(genotype=[100])
            ],
            preOps=sim.StepwiseMutator(rates=v0, incProb=incProb, loci=0),
            matingScheme=sim.RandomMating(),
            gen=t-t0
        )
        # make a copy of pop
        pop1 = pop.clone()
        # Evolve for another t0 generations with population size N1
        pop1.evolve(
            preOps=sim.StepwiseMutator(rates=v1, incProb=incProb, loci=0),
            matingScheme=sim.RandomMating(subPopSize=N1),
            gen=t0
        )
        # draw a random allele from pop1
        ind1 = pop1.individual(randint(0, N1-1))
        # if the allele length is greater than the threshold,
        if ind1.allele(randint(0,1)) > thresh:
            # Evolve pop2 for t0 generations with population size N2
            pop.evolve(
                preOps=sim.StepwiseMutator(rates=v2, incProb=incProb, loci=0),
                matingScheme=sim.RandomMating(subPopSize=N2),
                gen=t0
            )
            return (sum(pop1.genotype()) - sum(pop.genotype()))/(2.*N1)

