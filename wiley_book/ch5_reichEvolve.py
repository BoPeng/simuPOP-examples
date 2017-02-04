import simuOpt
simuOpt.setOptions(quiet=True, alleleType='long')
import simuPOP as sim
from reichDemo import demoModel
from reichStat import Ne
from simuPOP.utils import migrIslandRates
from itertools import product

def evolvePop(model, N0, N1, G0, G1, initSpec, mu, k, fitness,
        m=1, migrRate=0, logfile='', sp_logfile='', **kwargs):
    '''Evolve a population with specified allele frequencies (parameter
    initSpec) using given demographic (model, N0, N1, G0, G1, m), mutation
    (a k-allele model with parameters mu and k) and natural selection models
    (a multi-locus selection model with fitness vector s). Total disease
    allele frequency and effective number of alleles in the population
    and in all subpopulations are recorded if names of log files are provided.
    This function returns a tuple of these two statistics at the end of the
    evolution. Additional keyword arguments could be used to control when and
    how often statisitcs are outputed.
    '''
    L = len(fitness) // 3
    if not hasattr(mu, '__iter__'): # if a single mutation rate is given
        mu = [mu]*L
    # Create expressions to output f_e and ne at all loci, which are
    #   "%d\t%.4f\t%.4f\n" % (gen, 1-alleleFreq[x][0], ne[x])
    # for locus x.
    statExpr = '"%d' + r'\t%.4f\t%.4f'*L + r'\n" % (gen,' + \
        ', '.join(['1-alleleFreq[%d][0], ne[%d]' % (x, x) for x in range(L)]) + ')'
    demo_func = demoModel(model, N0, N1, G0, G1, m)
    pop = sim.Population(size=demo_func(0), loci=[1]*L,
            infoFields=['fitness', 'migrate_to'])
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitGenotype(freq=initSpec)
        ],
        preOps=[
            sim.KAlleleMutator(k=k, rates=mu, loci=range(L)),
            sim.MlSelector([
                sim.MaSelector(loci=i, fitness=fitness[3*i:3*(i+1)]) 
                for i in range(L)], mode=sim.MULTIPLICATIVE),
            sim.Migrator(rate=migrIslandRates(migrRate, m), begin=G0+1),
        ],
        matingScheme=sim.RandomMating(subPopSize=demo_func),
        postOps=[
            sim.IfElse(logfile != '' or sp_logfile != '', 
                Ne(loci=sim.ALL_AVAIL, vars=['ne'] if m == 1 else ['ne', 'ne_sp']),
                **kwargs),
            sim.IfElse(logfile != '', 
                sim.PyEval(statExpr, output='>>' + logfile), **kwargs),
            sim.IfElse(m > 1 and sp_logfile != '', 
                sim.PyEval(statExpr, output='>>' + sp_logfile,
                # subPops=sim.ALL_AVAIL will evalulate the expression in each
                # subpopulation's local namespace (vars(sp)).
                subPops=sim.ALL_AVAIL, begin=G0), **kwargs),
            ],
        finalOps=Ne(loci=sim.ALL_AVAIL),
        gen = G0 + G1
    )
    return tuple([1-pop.dvars().alleleFreq[x][0] for x in range(L)] + \
        [pop.dvars().ne[x] for x in range(L)])

