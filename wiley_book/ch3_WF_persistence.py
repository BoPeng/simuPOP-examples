import simuPOP as sim
from math import log
N = 100
p = 0.4
pop = sim.Population(2*N, ploidy=1, loci=1)
# Use a simulator to simulate 500 populations simultaneously.
simu = sim.Simulator(pop, rep=500)
gens = simu.evolve(
    initOps=sim.InitGenotype(prop=[p, 1-p]),
    # A RandomSelection mating scheme choose parents randomly regardless
    # of sex and copy parental genotype to offspring directly.
    matingScheme=sim.RandomSelection(),
    postOps=[
        # calculate allele frequency at locus 0
        sim.Stat(alleleFreq=0),
        # and terminate the evolution of a population if it has no
        # allele 0 or 1 at locus 0.
        sim.TerminateIf('alleleFreq[0][0] in (0, 1)'),
    ],
)
# find out populations with or without allele 1
gen_lost = []
gen_fixed = []
for gen,pop in zip(gens, simu.populations()):
    if pop.dvars().alleleFreq[0][0] == 0:
        gen_lost.append(gen)
    else:
        gen_fixed.append(gen)

print('''\nMean persistence time: %.2f (expected: %.2f)
Lost pops: %d (expected: %.1f), Mean persistence time: %.2f (expected: %.2f)
Fixed pops: %d (expected: %.1f), Mean persistence time: %.2f (expected: %.2f)'''\
% (float(sum(gens)) / len(gens), -4*N*((1-p)*log(1-p) + p*log(p)),
len(gen_fixed), 500*p, float(sum(gen_fixed)) / len(gen_fixed),
-4*N*(1-p)/p*log(1-p), len(gen_lost), 500*(1-p),
float(sum(gen_lost)) / len(gen_lost), -4*N*p/(1-p)*log(p)))

