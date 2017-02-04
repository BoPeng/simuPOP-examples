import simuOpt
simuOpt.setOptions(gui=False, alleleType='binary')
import simuPOP as sim
import random, math

from ch6_genCaseCtrl import genCaseControlSample

def penetrance(alpha, beta1, beta2, beta3, gamma1, gamma2):
    def func(geno):
        e = random.randint(0, 1)
        g1 = geno[0] + geno[1]
        g2 = geno[2] + geno[3]
        logit = alpha + beta1*g1 + beta2*g2 + beta3*g1*g2 + gamma1*e*g1 + gamma2*e*g2
        return 1 / (1. + math.exp(-logit))
    return func

sample = genCaseControlSample(pop, 1000, 1000,
    sim.PyPenetrance(func=penetrance(-5, 0.20, 0.4, 0.4, 0.2, 0.4),
        loci=['rs4491689', 'rs6869003']))
sim.stat(sample, association=sim.ALL_AVAIL)
# get p-values
sample.dvars().Allele_ChiSq_p
