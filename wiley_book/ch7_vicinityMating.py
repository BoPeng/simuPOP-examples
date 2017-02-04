import simuPOP as sim
from random import randint, uniform, normalvariate
class VicinityMating(sim.HomoMating):
    '''A homogeneous mating scheme that choose parents that are close to
    each other according to values at an information field.'''
    def __init__(self, locationField='x', varOfLocation=1, vicinity=1,
        numOffspring=1, sexMode=sim.RANDOM_SEX, ops=sim.MendelianGenoTransmitter(),
        subPopSize=[], subPops=sim.ALL_AVAIL, weight=0):
        '''Creates a random mating scheme that selects a parent randomly and
        another random parent who is in vivinity with him/her, namely with
        location that is within [x-v, x+v] where x is the location of the first
        parent, and v is specified by parameter vicinity. For each offspring,
        its location is set according to a normal distribution with a mean that
        is the average of parental locations, and a variance varOfLocation.
        '''
        self.field = locationField
        self.vicinity = vicinity
        self.varOfLocation = varOfLocation
        if hasattr(ops, '__iter__'): # if a sequence is given
            # WithArgs is needed because field name is a variable.
            allOps = ops + [sim.PyTagger(sim.WithArgs(self._passLocation, [self.field]))]
        else:
            allOps = [ops, sim.PyTagger(sim.WithArgs(self._passLocation, [self.field]))]            
        sim.HomoMating.__init__(self,
            chooser = sim.PyParentsChooser(self._chooseParents),
            generator = sim.OffspringGenerator(allOps, numOffspring, sexMode),
            subPopSize = subPopSize,
            subPops = subPops,
            weight = weight)
    
    def _passLocation(self, field):
        return normalvariate((field[0]+field[1])/2, self.varOfLocation)
    
    def _chooseParents(self, pop, subPop):
        # sort individuals according to location
        pop.sortIndividuals(self.field)
        while True:
            # select the first parent
            p1 = randint(0, pop.subPopSize(subPop) - 1)
            x1 = pop.individual(p1).info(self.field)
            s1 = pop.individual(p1).sex()
            # find all inviduals with opposite sex within vivinity of p1
            inds = []
            p = p1 + 1
            while p < pop.subPopSize(subPop) and \
                pop.individual(p, subPop).info(self.field) < x1 + self.vicinity:
                if pop.individual(p, subPop).sex() != s1:
                    inds.append(p)
                p += 1
            p = p1 - 1
            while p >= 0 and \
                pop.individual(p, subPop).info(self.field) > x1 - self.vicinity:
                if pop.individual(p, subPop).sex() != s1:
                    inds.append(p)
                p -= 1
            # if no one is invicinity, find another pair
            if len(inds) == 0:
                continue
            # choose another parent
            p2 = inds[randint(0, len(inds) -1)]
            # return indexes of both parents
            if s1 == sim.MALE:
                yield p1, p2
            else:
                yield p2, p1

pop = sim.Population(size=2000, loci=1, infoFields='x')
# define VSPs x<1, 1<=x<2, 2<=x<3, 3<=x<4, ...
pop.setVirtualSplitter(sim.InfoSplitter(field='x', cutoff=range(1, 8)))
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitInfo(lambda : uniform(0, 8), infoFields='x'),
        # only individuals in the middle range has certain genotype
        sim.InitGenotype(freq=[0.6, 0.4], subPops=[(0, 4)]),
    ],
    matingScheme=VicinityMating(locationField='x', vicinity=1, varOfLocation=0.5),
    postOps=[
        sim.Stat(alleleFreq=0, subPops=[(0, sim.ALL_AVAIL)], vars='alleleFreq_sp'),
        sim.PyEval(r"'%.3f ' % alleleFreq[0][1]", subPops=[(0, sim.ALL_AVAIL)]),
        sim.PyOutput('\n'),
    ],
    gen = 10
)    

