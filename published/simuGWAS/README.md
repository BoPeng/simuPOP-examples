(:title Simulation of samples for Genome-Wide Association Studies:)
%rfloat text-align=center margin-top=5px margin-right=25px margin-bottom=15px margin-left=25px % [[Attach:simuGWAS.zip | http://simupop.sourceforge.net/images/download.jpg]]|simuGWAS.zip

## Introduction

This script simulate GWAS data using data from the HapMap project. Please refer to 

->Peng B, Amos CI. (2010) [[http://www.biomedcentral.com/1471-2105/11/442|'''Forward-time simulation of realistic samples for genome-wide association studies''']] BMC Bioinformatics. 2010 Sep 1;11:442. doi: 10.1186/1471-2105-11-442.

and cite this article if you have used this script for your research.

## List of scripts
The file simuGWAS.zip contains several files, including:

:'''[[Attach:simuGWAS.py | simuGWAS.py]]''': This script evolves a population forward in time, subject to rapid population expansion, mutation, recombination and natural selection. A trajectory simulation method is used to control the allele frequency of optional disease predisposing loci. A scaling approach can be used to improve efficiency when weak, additive genetic factors are used

:'''[[Attach:loadHapMap2.py | loadHapMap2.py]]''': This python module provides function loadHapMapPop to download and import the HapMap populations. It can also be served as a script to download part or all populations from the phase 2 of the HapMap dataset

:'''[[Attach:loadHapMap3.py | loadHapMap3.py]]''': This script downloads and loads release 3 of hapmap phase 3 datasets in ftp://ftp.ncbi.nlm.nih.gov/hapmap//phasing/2009-02_phaseIII/HapMap3_r2/ and downloads the fine-scale recombination map from http://ftp.hapmap.org/recombination/2008-03_rel22_B36/rates/ and saves the genetic distance of each marker in a dictionary (geneticMap) in each population's local namespace.

->The saved populations have the following features:

->1. Different populations are saved in different files. These populations may not be merged directly because they have different set of markers. Subpopulation name is specified ('ASW', 'CEU', 'CHB'...).

->2. Chromosome names are saved as "1", "2", "3", ...

->3. Basepairs are used to specify physical distances of loci.

->4. Alleles are saved as 0 and 1 as appear in the HapMap datafile. Allele names such as 'A', 'G' are saved for each marker.

->5. A dictionary 'geneticMap' is used to store genetic distance of each marker.

:'''[[Attach:selectMarkers.py | selectMarkers.py]]''': This python module provides several utility functions that handles HapMap populations. When used as a script, this module creates a population using selected markers and populations

:'''[[Attach:example1.py | example1.py]]''': This example generates a sample using default parameters. and compares it with the HapMap population from which it is generated

:'''[[Attach:example2.py | example2.py]]''': This example simulates a case-control sample using a gene environment interaction model. 

->%red%note:%% This example requires a marker list file from Illumina. If you do not have this file, you can get a similar file from http://www.openbioinformatics.org/gengen/gengen_download.html. You should then change line 46 of example 2 from `ann = open('HumanHap550v3_A.lst')` to `ann = open('hh550v3_snptable.txt')`, and line 49 from `names.append(line.split(',')[1])` to `names.append(line.split('\t')[0])`.

:'''[[Attach:example3.py | example3.py]]''': This example simulates recent and remote selective sweep and draw trio samples

:'''[[Attach:example4.py | example4.py]]''': This example simulates an admixed population

 %color=red%These examples are described in detailed in [[http://www.biomedcentral.com/1471-2105/11/442/abstract|this paper]]%%. Please do not hesitate to contact me for any question. More examples will be added when simuGWAS.py is used to produce other types of samples.

## Steps of simulation

# Download HapMap data using loadHapMap2.py or loadHapMap3.py. Because of the size of the datasets, you do not have to download all the data at once.

# Select individuals and markers from HapMap 2 or 3 markers using selectMarkers.py. You can select markers by chromosome, region, number of markers, marker distance, and minor allele frequency, or from a list of markers stored in a file.

# Run simuGWAS.py to evolve the initial population, with optional disease loci and selection pressure. This step will generate a population with a large number of individuals.

# Following one of the examples to process the simulated population and draw samples. If your simulation closely mimics what one of the examples does, you can use the script with perhaps changes to the parameters. Otherwise, you will need to learn some simuPOP in order to change these script or write your own post-processing scripts.

## How to execute these scripts

!!!! Graphical user interface
The easiest way to run loadHapMap2.py, loadHapMap3.py, selectMarkers.py and simuGWAS.py is to use their graphical user interfaces. When you execute these scripts without parameter, parameter input dialogs will be displayed so that you can input parameters interactively. For example, the GUI of script simuGWAS.py is

Attach:simuGWAS.jpg

!!!! Using command line parameters

If you need to run these scripts in batch mode, you can use parameter `--gui=False` to disable GUI and use command line options to input parameters. If some of the parameters are not specified, your will be prompted to input them. Please use 

->`> simuGWAS.py -h`

to check what options are available.

!!!! Use a Python script.

All these scripts could be imported from another script. Because it is clear to specify all parameters in a script, all examples (example1.py, ... example4.py) import and executes script in this way. After these populations are simulated, these scripts use different techniques to analyze them.



!!!!%block class=messagehead%13 June 2011

>>message<<
!20:44 by '''[[~Matt Johnson]]'''!
It should be noted that to run the script example2.py as provided, the file HumanHap550v3_A.lst (list of all Illumina markers) is needed. I downloaded what I think is a similar file from here: http://www.openbioinformatics.org/gengen/gengen_download.html

I had to change the file to be comma delimited rather than tabs. This list of Illumina markers seems to be different, because 'rs4491689' is the 2944th marker on chromosome 2, and the original example2.py script only grabs 2000 markers on chromosome 2. Once I updated this to grab 3000 markers from each chromosome, example2.py executed correctly.

I can upload a corrected version if it's desired.
>><<


!!!!%block class=messagehead%09 August 2013

>>message<<
!11:18 by '''[[~Leslie Foldager]]'''!
Note that in simuGWAS.py lines 227-229 it is indicated that an exponential expansion model is used. But the paper (BMC Bioinformatics 2010) says that linear expansion was used and indeed the script produces a linear growth. The function is also identical to the 'linearExpansion' function in the source code 5.6 on page 125 in the book by Peng, Kimmel and Amos (ISBN: 9780470503485).

The reason for choosing 50000 as the default final population size is a bit unclear to me. In simuGWAS.py it is noted that this "is the recommended value when all HapMap populations are used (60+60+90)*200". Yes there are 60+60+90 unrelated subjects in HapMap2 but the examples in the paper use HapMap3 (993 unrelated individuals) and furthermore (60+60+90)*200=42000 ... so why 50000? Choosing 50000 actually results in an expected effective population size of 12657 (according to the output that I got from running the script) which is probably a reasonable size for a present day population and very close to the number 12658 mentioned in the paper and book. So the 10^5 noted in the paper and book (the book actually says 105 but that must be a typo) should probably have been 50000!?

So I'm wondering: is 50000 the recommended final size when using the 993 HapMap3 or when using the 210 HapMap2 or in both cases?
>><<



!!!!%block class=messagehead%04 April 2014

>>message<<
!04:49 by '''Bo Peng'''!
The expDemoFunc in simuGWAS.py is actually a typo. Both exponential and linear population growth models were tested. The exponential growth model was removed but the function name was not corrected. We apologize for the confusion. Note that simuPOP 1.1.2 provides ways to specify demographic models easily so a newer version of simuGWAS.py will remove these function altogether. 

Yes, 50,000 is recommended because it yields an effective population size of about 12,000. If you are using HapMap2, the effective size of the final population will be much smaller because of the smaller initial population size. This means HapMap2 simulated data can only be used to simulate smaller, probably isolated populations.
>><<

