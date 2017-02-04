import simuPOP as sim
from random import uniform, randint
from math import exp

class LC_model:
    def __init__(self, LC_beta0_male, LC_beta0_female, LC_a0, COPD_beta0,
        COPD_a0, G1_smoking_rate, rr_G2, rr_G3, rr_G4, rr_random,
            rr_LC_maleSmoker, rr_LC_femaleSmoker, rr_COPD_smoker,
            rr_LC_COPD):
        '''LC model with parameters for different relative risks
        A disease model is responsible for updating age_LC, age_death, and
        smoking. It can use other information fields to its own use.
        '''
        self.LC_beta0_male, self.LC_beta0_female, self.LC_a0, self.COPD_beta0, \
        self.COPD_a0, self.G1_smoking_rate, self.rr_G2, self.rr_G3, self.rr_G4, \
        self.rr_random, self.rr_LC_maleSmoker, self.rr_LC_femaleSmoker,  \
        self.rr_COPD_smoker, self.rr_LC_COPD = \
            LC_beta0_male, LC_beta0_female, LC_a0, COPD_beta0, \
            COPD_a0, G1_smoking_rate, rr_G2, rr_G3, rr_G4, rr_random, \
            rr_LC_maleSmoker, rr_LC_femaleSmoker, rr_COPD_smoker, rr_LC_COPD
           
  
    def _cdf(self, t, beta, a):
        # use a stepwise function to approximate integration
        ch = sum([beta*exp((i-18)*a) for i in range(18, t+1)])
        v = 1 - exp(-ch)
        return v
    
    def _ageOfOnset(self, beta, a):
        '''Calculate age of onset with a given beta and a. This is only used to initialize
        a population without considering cessation and other factors.
        '''
        u = uniform(0, 1)
        # a bisection method will be more efficient..
        for age in range(80):
            if self._cdf(age, beta, a) > u:
                return age
        return 100
     
    def initialize(self, ind):
        '''Determines the risk of LC for passed offividual. The return value is the
        multiple of a base hazard function.
        '''
        geno = [ind.allele(x,0) + ind.allele(x,1) for x in range(4)]
        # smoking, determined by G1
        ind.smoking = uniform(0,1) < self.G1_smoking_rate[geno[0]]
        # original age of death
        ind.age_death = randint(60, 80)
        # coefficient for LC and COPD
        LC_beta = self.LC_beta0_male if ind.sex() == sim.MALE else self.LC_beta0_female
        COPD_beta = self.COPD_beta0
        # smoking
        if ind.smoking:
            if ind.sex() == sim.MALE:
                LC_beta *= self.rr_LC_maleSmoker
            else:
                LC_beta *= self.rr_LC_femaleSmoker
            COPD_beta *= self.rr_COPD_smoker
        # G2
        LC_beta *= self.rr_G2[geno[1]]
        # G3
        if ind.smoking:
            LC_beta *= self.rr_G3[geno[2]]
        # G4
        COPD_beta *= self.rr_G4[geno[3]]
        if self._ageOfOnset(COPD_beta, self.COPD_a0) < ind.age_death:
            LC_beta *= self.rr_LC_COPD
        # random factor
        LC_beta *= 1 + self.rr_random * uniform(0,1)
        # LC?
        ind.age_LC = self._ageOfOnset(LC_beta, self.LC_a0)
        # adjust age of death if someone will get LC
        if ind.age_death < ind.age_LC + 6:
            ind.age_death = ind.age_LC + 6
        return ind.age <= ind.age_death
    
    def updateStatus(self, pop):
        # required by the evolutionary process but this disease model
        # currently does not need to update individual status dynamically.
        return True


