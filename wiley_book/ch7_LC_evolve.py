import simuPOP as sim

from ch7_LC_model import LC_model

def LC_evolve(popSize, alleleFreq, diseaseModel):
    '''
    '''
    pop = sim.Population(size=popSize, loci=[1]*len(alleleFreq),
        infoFields = ['age', 'smoking', 'age_death', 'age_LC', 'LC'])
    pop.setVirtualSplitter(sim.CombinedSplitter(splitters=[
        sim.InfoSplitter(field='age', cutoff=[20, 40],
            names=['youngster', 'adult', 'senior']),
        sim.SexSplitter(),
        sim.InfoSplitter(field='smoking', values=[0, 1, 2],
            names=['nonSmoker', 'smoker', 'formerSmoker'])
        ]
    ))
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitInfo(range(75), infoFields='age')] +
            [sim.InitGenotype(freq=[1-f, f], loci=i) for i,f in enumerate(alleleFreq)] + [
            sim.PyOperator(func=diseaseModel.initialize),
        ],
        preOps=[
            sim.InfoExec('age += 1'),
            # die of lung cancer or natural death
            sim.DiscardIf('age > age_death')
        ],
        matingScheme=sim.HeteroMating([
            sim.CloneMating(weight=-1),
            sim.RandomMating(ops = [
                sim.MendelianGenoTransmitter(), 
                sim.PyOperator(func=diseaseModel.initialize)],
                subPops=[(0, 'adult')])
            ],
            subPopSize=lambda pop: pop.popSize() + popSize/75),
        postOps = [
            # update individual, currently ding nothing.
            sim.PyOperator(func=diseaseModel.updateStatus),
            # determine if someone has LC at his or her age
            sim.InfoExec('LC = age >= age_LC'),
            # get statistics about COPD and LC prevalence
            sim.Stat(pop, meanOfInfo='LC', subPops=[(0, sim.ALL_AVAIL)],
                vars=['meanOfInfo', 'meanOfInfo_sp']),
            sim.PyEval(r"'Year %d: Overall %.2f%% M: %.2f%% F: %.2f%% "
                r"NS: %.1f%%, S: %.2f%%\n' % (gen, meanOfInfo['LC']*100, "
                r"subPop[(0,3)]['meanOfInfo']['LC']*100,"
                r"subPop[(0,4)]['meanOfInfo']['LC']*100,"
                r"subPop[(0,5)]['meanOfInfo']['LC']*100,"
                r"subPop[(0,6)]['meanOfInfo']['LC']*100)"),
        ],
        gen = 100
    )


if __name__ == '__main__':
    LC_evolve(10000, [0.5, 0.1, 0.2, 0.3], LC_model(
        LC_beta0_male=0.0025, LC_beta0_female=0.0015, LC_a0=0.012,
        COPD_beta0=0.00015, COPD_a0=0.01, G1_smoking_rate=[0.3, 0.4, 0.5],
        rr_G2=[1, 1.5, 1.8], rr_G3=[1, 1.1, 1.3], rr_G4=[1, 1.5, 2],
        rr_random=1.3, rr_LC_maleSmoker=10, rr_LC_femaleSmoker=8,
        rr_COPD_smoker=6, rr_LC_COPD=2))


