# Simulation of genome sequences with rare variants

## NEWS

* (Version 1.3) Jan, 13, 2014: The script has been updated to support simulation of regions on X-chromosome. A bug that prevents the use of some of the gamma distributions has been fixed.

* Feb, 3, 2012: A new distribution called mixed_gamma1 is added to srv.py that mimics the distribution used in Kyrukov et al, 2009. The major difference between mixed_gamma1 and mixed_gamma is that mixed_gamma discard fitness values outside of a range, whereas mixed_gamma1 returns 0 or 1 for values for values falling out of the boundaries.

* A bug that can cause incorrect genotype for recombinants under some rare conditions has been reported and fixed in trunk. If you are running the simulation with recombination, you should use simuPOP version 1.0.7 or higher.

* May 4, 2011: An updated version is uploaded. This version adds a parameter `postHook` to the function so that a Python function can be called, for example to draw a sample, at the end of each generation. An example can be seen [here](Attach:evol.py).

## NOTE

* This script is used to simulate mutants at all nucleotide locus over a short genome region, and should not be used to simulate sequences of, for example, more than 1 million basepairs (unless you are simulating a small population of no more than 10,000 individuals).

* If you notice an error for memory allocation error, it is likely that you are using a 32bit operating system. Please consider running this script on a 64 OS with 64 bit of python and simuPOP, and with at least 4G of RAM.

## Introduction

This script simulates the introduction and evolution of genetic variants in one or more ''regions'' of chromosomes. These regions span roughly `10k` to `100k` basepair and can be considered as a gene. During evolution, mutants are introduced to the population and change the fitness of individuals who carry these mutants. **The most distinguishing feature of this script is that it allows multi-locus fitness schemes with random or locus-specific diploid single-locus selection models to newly arising mutants**. A multi-locus selection model is used to assign a fitness value to individuals according the mutants they carry.

Please cite

  Bo Peng, Xiaoming Liu (2011) [Simulating Sequences of the Human Genome with Rare Variants](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3164177/) Hum Hered. 70(4): 287–291. Published online 2011 January 6. doi: 10.1159/000323316 PMCID: PMC3164177

If you have used srv for your research.

### Genotype structure

We assume one or more regions of chromosomes. Mutation can happen at any nucleotide locus which causes change of fitness of individuals carrying these mutants. The regions should be specified as `'ch1:1..50000'`. A list of regions is acceptable.

### Demographic model

This script uses a multi-stage population expansion / bottleneck model with population structure. Assuming there are n stages, the demographic model can be specified by

```
N = [N0, N1, ..., N_n-1, N_n]
G = [G0, G1, ..., G_n-1]
splitTo = [p1, p2, ..., pm]
and splitAt
```
where `N` is the starting population size and size at the end of each stage, `G` is the number of generations at each stage, `p1`, `p2`, ..., `pm` are proportions of subpopulations (should sum to 1). Then,

* If $N_t \lt N_{t+1}$, an exponential population expansion model is used to expand population from size {$N_t$} to {$N_{t+1}$}.
* If {$N_t > N_{t+1}$}, an instant population reduction model is used to reduce population size instantly to {$N_{t+1}$}. 
* If {$m > 1$}, the population will be split into `m` subpopulations according to proportions p1, ..., pm, at geneartion `splitAt`.

The default demographic model consists of a long burn-in stage, a short bottleneck stage and a rapid population expansion stage. The burn-in stage of this demographic model evolves a relatively small population of 8100 individuals until it reaches a mutation selection equilibrium. After a short bottleneck stage of 7900 individuals, the population grows exponentially to a population of 900,000 individuals in 370 generations. This demographic model reflects a demographic model of the European population (Kryukov, et al., 2009) and should be specified as

```
N=[8100, 8100, 7900, 900000], G=[8000, 100, 370]
````

### Mutation model

This script supports two mutation models

* A default `finite-sites mutation model` where mutations can happen at any locus. If a mutant is mutated, it will be mutated back to a wildtype allele.

* A `infinite-sites mutation model` where mutations can happen only at loci without existing mutant. If a mutation hits a locus with existing locus, it will be relocated to another locus without existing mutant, or be ignored if no such locus can be identified.

This script can output a mutant file that dump all mutation events during evolution. This file has the format of

```
generation location individual_index type
```

where type is 0 for forward, 1 for backward, 2 for relocated and 3 for ignored mutations.

%NOTE: If you are using an infinite-sites model, please pay attention to the number of segregation sites (column 3 of the output) and see if the population has been saturated with new mutants (#segregation site == #sites).

### Selection model

A random distribution is used to assign selection coefficients to newly arising mutations. This program incorporates multiple sets of selection parameters estimated from human genome data using different demographic models (Boyko, et al., 2008; Eyre-Walker, et al., 2006; Kryukov, et al., 2009; Williamson, et al., 2005).  For example, using a mixed gamma distribution, a mutant can have a selection coefficient of zeros (neutral alleles) or a random number drawn from a gamma distribution ranging from -0.00001 to -0.1 (Kryukov, et al., 2009).

Single locus selection model is specified using a selection coefficient ''s'' and a dominance coefficient ''h'' (default to 0.5), which assigns fitness 1, 1-hs and 1-s for genotypes 00, 01/10 and 11 (0 is wildtype allele) respectively.

If there are multiple mutants, the overall fitness of an individual is determined by either a multiplicative, an additive, or an exponential model. 

* An multiplicative model where {$f = \prod f_i$}
* An additive model where {$f=1-\sum\left(1-f_{i}\right)$}
* An exponential model (default) where {$f=\exp\left(\sum\left(1-f_{i}\right)\right)$}

NOTE: Fixed sites are reverted to wildtype alleles to avoid problems such as the Muller's Ratchet. You can remove operator RevertFixedSites from the script if this is not what you need.

### Recombination

Because this script simulates small regions of chromosomes, recombination is usually ignored. However, if you would like to simulate longer regions, or simulate special cases such as unlinked loci (r=0.5), you can specify a recombination rate using parameter `recRate`, which specifies recombination rate per basepair per generation.

This script works best for either full recombination (r=0.5) or realistic LD (e.g. r < 1e-5), specifying large recombination rates (e.g. `r=0.1`) will results in low performance due to excessive number of recombinations.

### Migration

If a population is split into m subpopulations, a migration rate `migrRate` can be specified to migrate individuals between subpopulations using an island model.

### Output

The end result of this script include

* A file that saves population statistics if parameter `statFile` is specified. Otherwise, the statistics will be written to standard output.

* A map file that contains the mutant location, frequency, selection coefficient and dominance coefficient. 

* A mutant file that contains location of mutants of each individual.

* One or more simuPOP population that can be imported to simuPOP and analyzed.

%blue%Note that fixed mutants are not counted as segregation sites but are included in the map and mutant files.

## Installation

This script requires simuPOP 1.0.5 to execute. The installation steps are described in detail in the simuPOP website. If you are using windows, please
* Download and install Python 2.6 from http://www.python.org.
* Download and install simuPOP 1.0.5 or later from http://sourceforge.net/projects/simupop/files/
* Download [srv.py](Attach:srv.py) and execute.

## How to use this script


### Run from a command line in batch mode

If you need to run it in batch mode, you can use command line, using options such as

```
simuRareVariants.py --gui=batch --selDist=gamma3
```

Default values will be used for unspecified parameters. 

### Import simuRareVariants from another script

You can import this script from another Python script and call its functions directly. Please see scripts in the examples section for details.

### Complete list of options

This is the output of `simuRareVariants.py -h`

```
Simulating a population of sequences forward in time, subject to mutation,
natural selection and population expansion. Because most mutants are introduced
during the rapid population expansion, most of the alleles will be rare at the
end of the simulation. Samples simulated using this script can be used to study
genetic diseases caused by a large number of rare variants.

usage: srv.py [--opt[=arg]] ...

options:
  -h, --help
        Display this help message and exit.

  --config=ARG (default: None)
        Load parameters from a configuration file ARG.

  --optimized
        Run the script using an optimized simuPOP module.

  --gui=[batch|interactive|True|Tkinter|wxPython] (default: None)
        Run the script in batch, interactive or GUI mode.

  --regions=ARG  (default: ['chr1:1..50000'])
        A region (in basepair) means a piece of chromosome in which mutations
        can happen. A region should be expressed as chrXX:YYYY..ZZZZ where XX is
        chromosome number, YYYY is the starting position in basepair and ZZZZ is
        the ending position in basepair. The starting position should be at
        least one. If multiple regions are specified as a list of regions, they
        are assumed to be unlinked and will segregate independently even if they
        are on the same chromosome.

  --initPop=ARG  (default: '')
        Name of an initial population. If this file exists, it will be loaded
        and the evolution will start from this population, instead of a blank
        population.

  --N=ARG  (default: [8100, 8100, 7900, 900000])
        Assuming a n stage demographic model, this parameter specifies
        population sizes at the beginning of evolution and at the end of each
        stage. N_0,...,N_n. If N_i < N_i+1, an exponential population expansion
        model will be used to expand population from size N_i to N_i+1. If N_i <
        N_i+1, an instant population reduction will reduce population size to
        N_i+1. For example
          N=(1000,1000,100,1000)
        simulates a three stage demographic model where a population of constant
        size goes through a bottleneck of 100 indiviudals, and then expands
        exponentially to a size of 1000.

  --G=ARG  (default: [5000, 10, 370])
        Numbers of generations of each stage of a n stage demographic model.
        This parameter should have n elements, in comparison to n+1 elements for
        parameter N.

  --splitTo=ARG  (default: [1])
        This parameter, if specified, should be a list of proportions that add
        up to 1. The length of this list specifies the number of subpopulations
        to split.

  --splitAt=ARG  (default: 0)
        Split the population at specified generation according to specified
        proportions.

  --mutationModel=ARG  (default: 'finite_sites')
        Mutation model. The default mutation model is a finite-site model that
        allows mutations at any locus. If a mutant is mutated, it will be
        mutated to a wildtype allele. Alternatively, an infinite-sites model can
        be simulated where new mutants must happen at loci without existing
        mutant, unless no vacant loci is available (a warning message will be
        printed in that case).

  --mu=ARG  (default: 1.8e-08)
        Mutation rate

  --selModel=ARG  (default: 'multiplicative')
        Multi-locus selection model, namely how to obtain an overall individual
        fitness after obtaining fitness values at all loci. This script supports
        three models:
            multiplicative: prod (f_i) Product of individual fitness.
            additive: max(0, 1 - sum(1-f_i)) One minus the combined selection
        deficiencies.
            exponential: exp(sum(1-f_i)) Exponential of combined selection
        deficiencies.
        Note that f_i can be equal to or greater than zero, which represents
        neutral loci, or loci under positive selection.

  --selDist=ARG  (default: 'constant')
        Distribution of selection coefficient for new mutants. Each distribution
        specifies s (selection coefficient) and h (dominance coefficient,
        default to 0.5 for additivity) that assign fitness values 1, 1-hs and
        1-s for genotypes AA (wildtype), Aa and aa, respectively. Note that
        positive s is used for negative selection so negative s is needed to
        specify positive selection. Note that we use 2k in the default
        distribution of Gamma distributions because theoretical estimates of s
        is for each mutant with 1-2s as fitness value for genotype aa in our
        model. This script currently support the following distributions:
        * constant: A single selection coefficient that gives each mutant a
        constant value s. The default parameter for this model is 0.01, 0.5. You
        can set selCoef to 0 to simulate neutral cases or a negative value for
        positive selection.
        * gamma1: A basic gamma distribution assuming a constant population size
        model (Eyre-Walker et al, 2006). The default parameters for this model
        is Pr(s=x)=Gamma(0.23, 0.185*2), with h=0.5. A scaling parameter 0.185*2
        is used because s in our simulation accounts for 2s for Eyre-Walker et
        al.
        * gamma2: A gamma distribution assuming a two-epoch population size
        change model for African population (Boyko et al, 2008). The default
        parameters for this model is Pr(s=x)=Gamma(0.184, 0.160*2), with h=0.5.
        * gamma3: A gamma distribution (for s) assuming a complex bottleneck
        model for European population (Boyko et al, 2008). The default
        parameters for this model is Pr(s=x)=Gamma(0.206, 0.146*2) with h=0.5.
        * mixed_gamma: Parameter of this model should be a list of (a, p, k,
        theta, h) where a is the probability of having s=p (neutral or adptive
        sites), k, theta are the parameter of a gamma distribution Recomended
        parameter is (0.0186, 0.0001, 0.184, 0.160*2, 0.5) for
        P(s=0.0001)=0.0186 and P(s=x)=(1-0.0186)*Gamma(0.184,0.160*2).
        If you would like to define your own selection model, please define your
        own function and pass it to parameter selDist of function
        simuRareVariants in the script.

  --selCoef=ARG  (default: None)
        Selection coefficient with its meaning determined by parameter selDist.
        If None is given, the default parameter for the selected distribution
        will be used. For example, parameter (0.001, 0) for a constant model
        defines a recessive model with fixed s. Note that a parameter of (k,
        theta, h) is needed for gamma distributions.

  --recRate=ARG  (default: 0)
        Recombination rate per base pair. If r times loci distance if greater
        than 0.5, a rate of 0.5 will be used.

  --migrRate=ARG  (default: 0)
        Migration rate to migrate individuals between subpoulations after the
        population is split into several subpopulations. An island model is
        used.

  --extMutantFile=ARG  (default: '')
        If a population is given, mutants from this population will be added to
        the population at specified generation. Only loci that are within the
        specified regions will be inserted. This population will be resized to
        population size at addMutantsAt before it is merged to the simulated
        population. This population is usually prepared using selectMarkers.py,
        using HapMap populations loaded using scripts loadHapMap2.py and
        loadHapMap3.py. These scripts are available from the simuPOP cookbook.

  --addMutantsAt=ARG  (default: 0)
        Generation number at which mutants from an external population will be
        inserted to the evolving population.

  --steps=ARG  (default: [100])
        Calculate and output statistics at intervals of specified number of
        generations. A single number or a list of numbers for each stage can be
        specified. If left unspecified, statistics at the beginning of each
        stage will be printed.

  --statFile=ARG  (default: '')
        File to output statistics. Default to standard output.

  --popFile=ARG  (default: 'output.pop')
        Filename to which the evolving population will be saved in simuPOP
        format. The default value of this parameter is 'output.pop', which saves
        the population at the end of the evolution to file 'output.pop'.
        Optionally, one or more generation numbers can be provided, in which
        case, the filename should be specified as an expression. For example,
        parameter ('!"output_%d.pop" % gen', (5000, -1)) saves the evolving
        population at the end of generation 5000, and the last generation.
        Please check the simuPOP user's guide for the use of expression in
        operator savePopulation.

  --markerFile=ARG  (default: 'output.map')
        Filename to which the marker information, including marker name
        (reg+index), chromosome, location, allele frequency and selection
        coefficient are saved. Monomorphic markers are ignored.

  --mutantFile=ARG  (default: 'output.mut')
        Filename to which the mutants are outputed. The file will be saved in
        the format of
            ind_id mut1 mut2 ...
        where ind_id is the index of individual (1, 2, ...), mut1 and mut2 are
        locations of mutants. Haplotypes for different regions and homologous
        chromosomes are saved in different lines in the order of
          reg1_ploidy1
          reg1_ploidy2
          reg2_ploidy1
          reg2_ploidy2
          ...

  --genotypeFile=ARG  (default: '')
        Filename to which the genotypes of all individuals are saved. The file
        will be saved in the format of
                famid id fa mo sex aff loc1_a1 loc1_a2 loc2_a1 loc2_a2 ...
        where famid is 1, 2, 3, ... id is always 1, fa, mo is always 0. Wildtype
        and mutant alleles are denoted by 0 and 1 respectively. This option is
        turned off by default because this format is not efficient in storing a
        small number of mutants.

  --verbose=ARG  (default: 1)
        0 for quiet, 1 for regular output, 2 for debug output. In the debug
        output, a file 'mutations.lst' will be saved with all mutation events.
        This option is not visible from gui.
```


### Examples

### Apply a penetrance or quantitative trait model

The selected population can be imported and post-processed using simuPOP. If you need to apply a quantitative trait model to the simulated population, you can use a function `pyQuanTrait` using a user-defined function. The only difference is that 'alleles' in the simulated population are locations of mutants. Penetrance model can be assigned similarly. 

Example [quanTraits.py](Attach:quanTraits.py) demonstrates how to apply a quantitative trait model and draw samples with extreme trait values.

Example [pedigree.py](Attach:pedigree.py) demonstrates how to evolve the simulated population for three more generations and draw three-generational pedigrees from the simulated multi-generational population, with restrictions on pedigree size and number of affected members.

These scripts import that calls the simuRareVariant function directly and uses functions defined in that script to save samples.

### Arbitrary distribution of selection coefficients.

The graphical user interface allows you to perform simulations for fix types (constant, gamma and mixed gamma). If you would like to define you own distribution, you can define a function that returns selection coefficient according to your distribution and pass it to the simuRareVariants function of simuRareVariants.py. You might need to use simuPOP's random number generation functions listed in [ here](http://simupop.sourceforge.net/manual_release/build/refManual_ch2_sec5.html#class-rng ).

Example [myDist.py](Attach:myDist.py) demonstrates how to define such a function.

### Location-specific selection coefficients

If you would like to define a selection model with selection coefficients related to mutation location. You can add a parameter `loc` to the distribution function. In that case, the location (in basepair) of the new mutant will be passed to your function. This feature allows you to define a neutral region within a larger region under selection, or return neutral for mutation happens at the last nucleotide of a codon. Moreover, if you have fixed set of selection coefficients, you can use this feature to pass them to the script.

Example [locSpecific.py](Attach:locSpecific.py) demonstrates how to define a fitness function that returns location-specific fitness values.

### Analyze all mutation events

If you would like to have a list of all mutation events happened during the evolutionary process, you can add `--verbose=2` to the command line. This will generate a file named `mutations.lst` which lists all mutations in the format of

->`generation  location  individual_index type`

where type is 0 for forward mutation, 1 for backward mutation and 2 for ignored mutation in the infinite-site model. You can process this file to trace the age of all mutants.


Example [mutAge.py](Attach:mutAge.py) demonstrates how to process this file and calculate the age of all mutants.

## Extending simuRareVariants.py

### Output more statistics such as the fitness value of everyone.

The default output of this script includes population size, number of segregation sites, average number of segregation sites per individual, average frequency of all mutants, average and mean fitness of individuals. You can modify the script to get more output. For example, you can use operator `infoEval` with `output='>>fitness.txt'"` to output fitness values for all individuals to a file. Please refer to the simuPOP user's guide on how to use these parameters.

## References

* Boyko AR, Williamson SH, Indap AR, Degenhardt JD, Hernandez RD, Lohmueller KE, Adams MD, Schmidt S, Sninsky JJ, Sunyaev SR, White TJ, Nielsen R, Clark AG, Bustamante CD: Assessing the evolutionary impact of amino acid mutations in the human genome. PLoS Genet 2008;4:e1000083.

* Eyre-Walker A, Woolfit M, Phelps T: The distribution of fitness effects of new deleterious amino acid mutations in humans. Genetics 2006;173:891-900.

* Kryukov GV, Pennacchio LA, Sunyaev SR: Most rare missense alleles are deleterious in humans: Implications for complex disease and association studies. Am J Hum Genet 2007;80:727-739. 

* Williamson SH, Hernandez R, Fledel-Alon A, Zhu L, Nielsen R, Bustamante CD: Simultaneous inference of selection and population growth from patterns of variation in the human genome. Proc Natl Acad Sci U S A 2005;102:7882-7887.
