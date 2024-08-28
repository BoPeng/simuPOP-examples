#!/usr/bin/env python
#script used in 'Roux &al, Molecular Biology and Evolution (2012)'
import simuOpt, types, os, sys, time
simuOpt.setOptions(alleleType='long')
import simuPOP as sim
from simuPOP import *
from simuPOP.utils import *
from simuPOP.sampling import drawRandomSample

try:
	from simuPOP.plotter import VarPlotter
except:
	useRPy=False
else:
	useRPy=True

options = [
    {
     'name':'NA',
     'default':5000,
     'label':'Ancestral Population Size',
     'type':[int, int],
     'validator':simuOpt.valueGT(0),
     'description':'Ancestral population size'
    },
    {
     'name':'N1',
     'default':1000,
     'label':'Daugther 1 Population Size',
     'type':[int, int],
     'validator':simuOpt.valueGT(0),
     'description':'Daughter 1 population size'
    },
    {
     'name':'N2',
     'default':1000,
     'label':'Daugther 2 Population Size',
     'type':[int, int],
     'validator':simuOpt.valueGT(0),
     'description':'Daughter 2 population size'
    },
     {
     'name':'Tbeforesplit',
     'default':1000,
     'label':'Number of Generations to simulate before split',
     'type':[int, int],
     'validator':simuOpt.valueGT(0),
     'description':'Number of Generations to simulate before split'
    },
    {
     'name':'Taftersplit',
     'default':200,
     'type':[int, int],
     'label':'Number of Generations to simulate after the split',
     'description':'Number of Generations to simulate after the split',
     'validator':simuOpt.valueGT(0)
    },
    {
     'name':'r2loci',
     'default':0.01,
     'label':'Recombination Rate between selected and neutral loci',
     'type':[float],
     'description':'Recombination rate',
     'validator':simuOpt.valueBetween(0., 1.),
    },
    {
     'name':'numLoci',
     'default':2,
     'label':'Number of loci',
     'type':[int],
     'description':'Number of loci',
     'validator':simuOpt.valueGT(0),
    },
	{
     'name':'K_sel',
     'default':10,
     'label':'Number of allelic states at selected locus',
     'type':[int, int],
     'validator':simuOpt.valueGT(0),
     'description':'Number of allelic states at selected locus'
    },
    {
     'name':'Mu_sel',
     'default':0.000001,
     'label':'Mutation rate at selected locus',
     'type':[float],
     'description':'Mutation rate at selected locus',
     'validator':simuOpt.valueBetween(0., 1.),
    },
    {
     'name':'S_sel',
     'default':0.1,
     'label':'Selection coefficient against homozygotes at selected locus',
     'type':[float],
     'description':'Selection coefficient against homozygotes at selected locus',
     'validator':simuOpt.valueBetween(0., 1.),
    },
	{
     'name':'L_neut',
     'default':1000,
     'label':'Number of biallelic positions in the neutral locus',
     'type':[int, int],
     'validator':simuOpt.valueGT(0),
     'description':'Number of biallelic positions in the neutral locus'
    },
   {
     'name':'Mu_neut',
     'default':0.000001,
     'label':'Mutation rate at each site of neutral locus',
     'type':[float],
     'description':'Mutation rate at each site of neutral locus',
     'validator':simuOpt.valueBetween(0., 1.),
    },
	{
     'name':'Nsample1',
     'default':50,
     'label':'Number of individuals to sample from population 1',
     'type':[int, int],
     'validator':simuOpt.valueGT(0),
     'description':'Number of individuals to sample from population 1'
    },
	{
     'name':'Nsample2',
     'default':50,
     'label':'Number of individuals to sample from population 2',
     'type':[int, int],
     'validator':simuOpt.valueGT(0),
     'description':'Number of individuals to sample from population 2'
    }
 ]

def simulate(NA, N1, N2, Tbeforesplit, Taftersplit, r2loci, numLoci, K_sel, Mu_sel, S_sel, L_neut, Mu_neut, Nsample1, Nsample2):
	pop = Population(size=NA, ploidy=2, loci=[L_neut+1], infoFields='fitness')
	def getfitness(geno):
     #  returns fitness of genotype geno at the overdominant locus with selection coeff S_sel
     # geno is (A1 A2)
		if geno[0] == geno[1] :
			return 1 - S_sel  # homozygote
		else:
			return 1       # heterozygote
	if useRPy:
		plotter=VarPlotter('alleleFreq[0][0],alleleFreq[0][1],alleleFreq[0][2],alleleFreq[0][3],alleleFreq[0][4]', ylim=[0,1], ylab='allele frequency', update=Tbeforesplit+Taftersplit-1, saveAs='slocus.png')
	else:
		plotter=NoneOp()
	g = pop.evolve(
		initOps = [
			InitSex(),
		#initially put 5 alleles at the selected locus with equal frequencies
			InitGenotype(loci=0, freq=[.1] * 10)
#			InitGenotype(loci=0, freq=[0.01, 0.1, 0.4, 0.2, 0.29])
			],
	preOps = [
		#Resize the ancestral population at the time immediatly before the split
		sim.ResizeSubPops([0], sizes=[N1+N2], at=Tbeforesplit-1),
	    	# split ancestral population in 3 subpopulations only works if NA>N1+N2
		sim.SplitSubPops(subPops=0, sizes=[N1, N2], at=Tbeforesplit),
		# apply overdominant selection by invoking function getfitness
		PySelector(loci=0, func=getfitness),
		],
	matingScheme = RandomMating(ops= [
		#apply recombination between the selected locus and the neutral locus at rate r2loci
		Recombinator(rates=r2loci, loci=0),
		]),
	postOps = [
		# apply mutation to the selected locus according to K allele model
		KAlleleMutator(k=K_sel, rates=[Mu_sel], loci=[0]),
		# apply mutation to the neutral sequence
		SNPMutator(u=Mu_neut,v=0,loci=list(range(1,L_neut))),
		#Computes the frequency of each allele at selected locus at the last generation
		Stat(alleleFreq=0,step=1),
		#output to the screen the frequency of the 4 first alleles at the selected locus
		PyEval(r'"%.0f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\n" % (gen, alleleFreq[0][0], alleleFreq[0][1], alleleFreq[0][2], alleleFreq[0][3], alleleFreq[0][4])',step=100),
		plotter,
		],
	# sets the last generation = Tbeforesplit+Taftersplit
	gen = Tbeforesplit + Taftersplit
	)
	#draw two random samples from species 1 and 2 with size Nsample1 and Nsample2
	sample = drawRandomSample(pop, sizes=[Nsample1,Nsample2])
	#write to file the content of the two random samples
	sim.utils.saveCSV(sample, filename='output.txt')
	return g


if __name__ == '__main__':
	pars = simuOpt.Params(options, doc ='This script simulates a speciation with variation at a selected locus under overdominant selection',details = '__doc__')


if not pars.getParam():
	sys.exit(0)


simulate(pars.NA, pars.N1, pars.N2, pars.Tbeforesplit, pars.Taftersplit, pars.r2loci, pars.numLoci, pars.K_sel, pars.Mu_sel,pars.S_sel, pars.L_neut, pars.Mu_neut, pars.Nsample1, pars.Nsample2)
