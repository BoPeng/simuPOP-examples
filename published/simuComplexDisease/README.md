(:title Simulation of human genetic diseases with selection:)
%rfloat text-align=center margin-top=5px margin-right=25px margin-bottom=15px margin-left=25px % [[Attach:simuComplexDisease.py | http://simupop.sourceforge.net/images/download.jpg]]|simuComplexDisease.py


This script uses a significantly different mechanism to control the allele 
frequency of disease susceptibility loci than simuForward.py. I will describe 
the method briefly here. More details please see 

    Peng B, Amos CI, Kimmel M (2007) Forward-Time Simulations of Human 
    Populations with Complex Diseases. PLoS Genet 3(3): e47

This program simulates the evolution of a complex common disease under the 
influence of mutation, migration, recombination and population size change. 
Starting from a small founder population, each simulation will go through
the following steps:

* Simulate the trajectory of allele frequency using specified disease model.
* Burn-in the population with mutation and recombination
* Introduce disease alleles and evolve the population with pre-simulated  allele frequency.
* Population structure and migration are specified along with demographic models. The population can be split into several equally-sized subpopulations and then evolve independently, or with migration. 

The result of the simulation is a large multi-generation population. To analyze 
the population, you will typically need to 

* Apply a penetrance function to the population and determine the affectedness  for each individual
* Draw Population and/or pedigree based samples and save in popular formats so that the samples can be analyzed by other software like genehunter.

Another script analComplexDisease gives excellent examples on how to perform these tasks.
Please follow that script to perform your own analyses.

The program is written in Python using the simuPOP modules. For more information,
please visit simuPOP website http://simupop.sourceforge.net .


### Genotype structure and disease

With typical settings, each individual will have 10 chromosomes, each having
20 equal spaced microsatellite or SNP markers. A disease will be caused by 
several disease susceptibility loci (DSL) between markers. For example, a
DSL may be .3 unit to the right of marker 25 (the six marker on the second
chromosome). Since we assume that fitness is only determined by genotype, 
not affectedness status or trait value, we do not assign individual 
affectedness till the sampling stage.


### Evolutionary Scenario

The evolutionary process can be divided into three stages:

!!!! Burn-in stage

A founder population will be initialized with a small number of haplotypes.
All DSL will be initialized with wild type alleles ( no disease). This leads
to complete linkage disequilibrium between markers and no LD with DSL.
The population will evolve for a number of generations, subject to mutation 
and recombination. This will break down LD between markers and let (hopefully)
each marker reach a mutation-drift equilibrium.

The mutation happens only at non-DSL markers, following a symmetric
stepwise mutation model for microsattelite and a Juke-Cantor model for
SNP markers. Recombination is uniform across all markers.


!!!! No-migration stage


The population will be split into 10 subpopulations and starts to grow,
aiming at 100,000 or more individuals at the final generation. No migration 
is allowed between subpopulations so population heterogeneity will build up.

!!!! Mixing stage


Individuals from different subpopulations will be able to migrate following
a circular step-stone model. Population heterogeneity will be reduced to
a level depending on migration rate and length of this stage.

!!!!Introduction of disease


The disease allele frequency is simulated before the simulation is performed.
A single disease mutant is introduce to each DSL at simulated mutant-introduction
generation. The allele frequency then evolve according to the simulated frequency
trajectory.

However, you can also specify some free DSL, who will evolve in a different manner.
Namely, they will be brought to high allele frequency in your specified time frame,
and then evolve freely. This method may fail due to extinction of disease alleles,
but it has the advantage of being able to simulate linked DSL.


### Statistics Monitored

A number of statistics will be measured and saved. They are:

* Fst before and after mixing
* Observed heterogeneity before and after mixing
* Disease allele frequency trajectory.

You can add more stat() operator if you need to preserve more information.

