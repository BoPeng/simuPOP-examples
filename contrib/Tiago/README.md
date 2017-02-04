(:title Landscape genetics example:)
%rfloat text-align=center margin-top=5px margin-right=25px margin-bottom=15px margin-left=25px % [[Attach:land.py | http://simupop.sourceforge.net/images/download.jpg]]|land.py

The example attached provides a landscape genetics simulation, the main points are:

#2D map varying from coordinates 0 to 100 on each dimension

#Map imposes selection conditions (left is "cold", right is "hot". Bottom is "low", top is "high")

#Individuals have a set of loci that are adaptative to temperature and height.

#The population has a constant size of 1000

#Individuals are able to (random) mate in a certain radius

#Offspring are able to walk from the origin point (equidistant between father and mother). There is probability space of where to walk, making the individual more prone to walk in the direction to which it is adapted

#The initial population is adapted to cold, low (i.e., bottom-left corner). Random mutation introduces adaptation to new environments.

#For each variable (temperature, height) there are a set of loci connected. The more mutations on the set, the more able is the individual to move to hotter and higher zones)

#Mutations might be wild-recessive (ie, having a single mutant allele makes it prone to move away for the original corner), wild-dominant (needs both alleles to make an effect) or random (an allele is choosen at random to define the phenotype).

#Selection doesn't use simuPOP operators, instead there is joint probability of an individual being killed before mating (joint of the adaptation of being suited to the point in space where it is)

#Some simuPOP programming techniques might not be the best

!!Video

From the output of the script above you can generate a set of images with [[Attach:plot.py |this script]].

Here is a video example:

(:youtube Al_0BlSxcHI:)

The reason the individuals go from one corner to the opposite is two-fold:


#All loci are wild-recessive (as soon as you have a mutation you will move). This is probably realistic
#There is no notion of "crowding". All individuals can occupy the same point in space. This is less realistic, of course...

The color of individuals is revealing of the genotype: black is wild-type. More green means more "hot" (i.e., tends to go to the right). More blue means more "high" (i.e., tends to go up)



!!Comments on the code supplied

=python [=
pop = population(POPSIZE, 2,
    loci=[1]*LOCIPERDIM*2,
    infoFields = ['x', 'y']
)
=]

Note the infoFields x and y.


=python [=
def geoChooser(pop, sp):
    while (True):
        gender = 2
        while gender != 1:
            g1Idx = randint(0, pop.popSize()-1)
            gender = pop.individual(g1Idx).sex()
        g1 = pop.individual(g1Idx)
        g1x = g1.info('x')
        g1y = g1.info('y')
        g2Cases = []
        for idx in range(0, pop.popSize()-1):
            tempG2 = pop.individual(idx)
            if tempG2.sex()==2:
                tg2x = tempG2.info('x')
                tg2y = tempG2.info('y')
                dist =  sqrt((g1x-tg2x)**2 + (g1y-tg2y)**2)
                if dist<=MATERADIUS: g2Cases.append(idx)
        if len(g2Cases)==0: continue
        yield g1Idx, g2Cases[randint(0,len(g2Cases)-1)]
=]

This is the random mating (but only within a certain radius).

=python [=
def placeIndividual(pop, off, dad, mom):
    cX = (dad.info('x') + mom.info('x'))/2
    cY = (dad.info('y') + mom.info('y'))/2
    #This is really a SQUARE
    x=-1
    y=-1
    xMut, yMut = countMutants(off)
    optX = getOptimalCoord(maxX, xMut)
    optY = getOptimalCoord(maxY, yMut)
    while x<0 or x>maxX: # or x<cX-WALKRADIUS or x>cX+WALKRADIUS:
        #We will assure that the most distant extreme is within 68?%
        d1X = abs((cX-WALKRADIUS)-optX) 
        d2X = abs((cX+WALKRADIUS)-optX)
        if d1X>d2X:
            bigger=cX-WALKRADIUS
            smaller=cX+WALKRADIUS
        else:
            bigger=cX+WALKRADIUS
            smaller=cX-WALKRADIUS
        if smaller<0: smaller=0
        if bigger<0: bigger=maxX
        x = normalvariate(smaller, WALKRADIUS)
    while y<0 or y>maxY: # or y<cY-WALKRADIUS or y>cY+WALKRADIUS:
        #We will assure that the most distant extreme is within 68?%
        d1Y = abs((cY-WALKRADIUS)-optY) 
        d2Y = abs((cY+WALKRADIUS)-optY)
        if d1Y>d2Y:
            bigger=cY-WALKRADIUS
            smaller=cY+WALKRADIUS
        else:
            bigger=cY+WALKRADIUS
            smaller=cY-WALKRADIUS
        if smaller<0: smaller=0
        if bigger<0: bigger=maxY
        y = normalvariate(smaller, WALKRADIUS)

    off.setInfo(x, pop.infoIdx('x'))
    off.setInfo(y, pop.infoIdx('y'))
    return True

=]

This places the new individual as a function of parent position and approximation to the point in space to which it is most adapted.

=python [=

def countMutants(ind):
    xMut = 0
    yMut = 0
    for i in range(LOCIPERDIM):
        if epiState==0: #wild-recissive
            if (ind.allele(i) == 1 or ind.allele(i+LOCIPERDIM*2) == 1): xMut+=1
            if (ind.allele(i+LOCIPERDIM) == 1 or ind.allele(i+LOCIPERDIM+LOCIPERDIM*2) == 1):
                yMut+=1
        elif epiState==1: #wild-dominant
            if (ind.allele(i) == 1 and ind.allele(i+LOCIPERDIM*2) == 1): xMut+=1
            if (ind.allele(i+LOCIPERDIM) == 1 and ind.allele(i+LOCIPERDIM+LOCIPERDIM*2) == 1):
                yMut+=1
        elif epiState ==2: #random
            xMut += ind.allele(i+randint(0,1)*2*LOCIPERDIM)
            yMut += ind.allele(i+LOCIPERDIM+randint(0,1)*2*LOCIPERDIM)
    return xMut, yMut


def getPenalty(maxPos, pos, numMuts):
    return 1.0 - abs(pos - getOptimalCoord(maxPos, numMuts))/maxPos

def killUnfit(pop):
    indivs = []
    for indPos in range(pop.popSize()):
        ind     = pop.individual(indPos)
        x       = ind.info('x')
        y       = ind.info('y')
        xMut, yMut = countMutants(ind)
        xFactor = getPenalty(maxX, x, xMut) 
        yFactor = getPenalty(maxY, y, yMut) 
        luck    = uniform(0,1)
        if luck>xFactor*yFactor:
            indivs.append(indPos)
    pop.removeIndividuals(indivs)
    return True
=]

This kills an unfit individual. Each coordinate penalty is calculated as a distance from the optimal point (which is dependent on the number of mutations). Penalties are multiplied (one for each coordinate) and compared to a random number.

=python [=
pop.evolve(
    initOps = [
        InitGenotype(freq=[1,0])
    ],
    preOps = [
        PyOperator(killUnfit),
        KamMutator(k=2, rates=[0.01]*LOCIPERDIM*2, loci=range(LOCIPERDIM*2)),
    ],
    matingScheme = HomoMating(
        PyParentsChooser(geoChooser),
        OffspringGenerator(ops=[
            MendelianGenoTransmitter(),#numOffspring=(UniformDistribution, 2, 4)),
            PyOperator(placeIndividual)]),
        subPopSize=[POPSIZE]
    ),
    postOps = PyOperator(dmp),
    gen = 1000
)
=]

The main loop, note the fixed initial frequencies, the mutation, placing the new individuals and the killings.


!!Copyright

The author of this script is [[http://tiago.org | Tiago Antao]]. You are free to use this script, but I would appreciate if you inform me in case you use it. In case of doubts, ideas, suggestions or whatever, feel free to contact me.

!!!!%block class=messagehead%10 September 2010

>>message<<
!13:14 by '''anonymous'''!"KamMutator not defined..."

I am a beginner both in SimuPop and Python so may explain why ...
>><<

