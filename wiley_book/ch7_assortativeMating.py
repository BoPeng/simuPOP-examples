import simuPOP as sim
from random import normalvariate
sigma = 1
def traits(geno):
    'genotypes are arranged as a1a2b1b2c1c2... where a,b,c are specified loci'
    A = sum(geno[:20]) + normalvariate(0, 2.5)
    B = sum(geno[20:40]) + normalvariate(0, 2.5)
    I = sum(geno[40:60]) + normalvariate(0, 2.5)
    D = B + I - A + normalvariate(0, sigma**2)
    return A, B, I, D

pop = sim.Population(100000, loci=[1]*40, infoFields=['A', 'B', 'I', 'D'])
pop.evolve(
    initOps=[
        sim.InitSex(maleProp=0.5),
        sim.InitGenotype(freq=[0.5, 0.5]),
    ],
    preOps=[
        sim.PyQuanTrait(func=traits, loci=sim.ALL_AVAIL,
            infoFields=['A', 'B', 'I', 'D']),
        sim.PyOperator(func=lambda pop: pop.sortIndividuals('D') is None),
    ],
    matingScheme=sim.HomoMating(
        chooser=sim.SequentialParentsChooser(),
        generator=sim.OffspringGenerator(
            ops=sim.MendelianGenoTransmitter(),
            numOffspring=2, sexMode=(sim.NUM_OF_MALES, 1))
    ),
    finalOps=sim.PyQuanTrait(func=traits, loci=sim.ALL_AVAIL,
            infoFields=['A', 'B', 'I', 'D']),
    gen=10
)

from rpy import r
def genoTraitCorrelation(loc, trait):
    'Calculate correlation between trait and genotype at a locus'
    geno = [ind.allele(loc,0) + ind.allele(loc,1) for ind in pop.individuals()]
    qtrait = pop.indInfo(trait)
    return r.cor(geno, qtrait)

# correlation between genotype at A loci with trait A
AA = [genoTraitCorrelation(loc, 'A') for loc in range(10)]
print(', '.join(['%.3f' % abs(x) for x in AA]))
# correlation between genotype at A loci with trait B (spurious)
AB = [genoTraitCorrelation(loc, 'B') for loc in range(10)]
print(', '.join(['%.3f' % abs(x) for x in AB]))
# correlation between genotype at unrelated loci with trait A
UA = [genoTraitCorrelation(loc, 'A') for loc in range(30, 40)]
print(', '.join(['%.3f' % abs(x) for x in UA]))
