import simuPOP as sim
from simuPOP.utils import simulateForwardTrajectory, simulateBackwardTrajectory, \
    migrSteppingStoneRates

def linearExpansion(N0, N1, G):
    '''Return a linear population expansion demographic function that expands
    a population from size N0 to N1 linearly in G generations. N0 and N1 should
    be a list of subpopulation sizes.'''
    step = [float(x-y) / G for x,y in zip(N1, N0)]
    def func(gen):
        if gen == G - 1:
            return N1
        return [int(x + (gen + 1) * y) for x, y in zip(N0, step)]
    return func

def simuGWAS(pop, mutaRate=1.8e-8, recIntensity=1e-8, migrRate=0.0001,
    expandGen=500, expandSize=[10000], DPL=[], curFreq=[], fitness=[1,1,1],
    scale=1, logger=None):
    # handling scaling...
    mutaRate *= scale
    recIntensity *= scale
    migrRate *= scale
    expandGen = int(expandGen / scale)
    fitness = [1 + (x-1) * scale for x in fitness]
    pop.dvars().scale = scale
    # Demographic function
    demoFunc = linearExpansion(pop.subPopSizes(), expandSize, expandGen)
    # define a trajectory function
    trajFunc = None
    introOps = []
    if len(DPL) > 0:
        stat(pop, alleleFreq=DPL, vars='alleleFreq_sp')
        currentFreq = []
        for sp in range(pop.numSubPop()):
            for loc in pop.lociByNames(DPL):
                currentFreq.append(pop.dvars(sp).alleleFreq[loc][1])
        # if there is no existing mutants at DPL
        if sum(currentFreq) == 0.:
            endFreq=[(x-min(0.01,x/5.), x+min(0.01, x/5., (1-x)/5.)) for x in curFreq]
            traj=simulateForwardTrajectory(N=demoFunc, beginGen=0, endGen=expandGen,
                beginFreq=currentFreq, endFreq=endFreq, nLoci=len(DPL),
                fitness=fitness, maxAttempts=1000, logger=logger)
            introOps=[]
        else:
            traj=simulateBackwardTrajectory(N=demoFunc, endGen=expandGen, endFreq=curFreq,
                nLoci=len(DPL), fitness=fitness, minMutAge=1, maxMutAge=expandGen,
                logger=logger)
            introOps = traj.mutators(loci=DPL)
        if traj is None:
            raise SystemError('Failed to generated trajectory after 1000 attempts.')
        trajFunc=traj.func()
    if pop.numSubPop() > 1:
        pop.addInfoFields('migrate_to')
    pop.dvars().scale = scale
    pop.evolve(
        initOps=sim.InitSex(),
        preOps=[
            sim.SNPMutator(u=mutaRate, v=mutaRate),
            sim.IfElse(pop.numSubPop() > 1,
                sim.Migrator(rate=migrSteppingStoneRates(migrRate, pop.numSubPop()))),
            ] + introOps,
        matingScheme=sim.ControlledRandomMating(loci=DPL, alleles=[1]*len(DPL),
            freqFunc=trajFunc, ops=sim.Recombinator(intensity=recIntensity),
            subPopSize=demoFunc),
        postOps = [
            sim.Stat(popSize = True, structure=range(pop.totNumLoci())),
            sim.PyEval(r'"After %3d generations, size=%s\n" % ((gen + 1 )* scale, subPopSize)'),
            sim.IfElse(pop.numSubPop() > 1,
                sim.PyEval(r"'F_st = %.3f\n' % F_st", step=10), step=10),
        ],
        gen = expandGen
    )
    return pop

