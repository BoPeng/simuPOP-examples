#  A power calculator for case control association studies of samples with  known family histories

%red%Due to a regression bug in simuPOP 1.0.4 and 1.0.5, please use simuPOP 1.0.3 or simuPOP 1.0.6 or later to execute this script.%%



### ChangeLog

* Jan, 30, 2012: Modify the 'Additive_AA' model to use an alternative definition of detectable relative risk. Use a slightly different power analysis formula.


### Introduction
This program is a power calculator for case-control association studies. Compared
to other power calculators, this program is unique in that it assumes that we
know the family history of cases and/or controls. Because cases with known family
history have higher probability of carrying a disease predisposing loci, such
study design is more powerful than regular case control association studies.

Please cite the following paper if you have used this script for your analysis:

->'''Power analysis for case–control association studies of samples with known family histories'''
->Bo Peng, Biao Li, Younghun Han and Christopher I. Amos
->LINK: http://www.springerlink.com/content/y85k314765524153/

### Installation

Because this script is written in Python and uses simuPOP for its graphical user interface, you will need to install Python (2.4 or higher, but not Python 3) and [ simuPOP](http://simupop.sf.net ) in order to use it. More specifically, you need to 

* Install Python if it is not available.
* Download and install %red%simuPOP 1.0.3%% or 1.0.6 or later following the installation instruction [ here](http://simupop.sourceforge.net/Main/Download ).
* Download caseControlPower.py from this page and save it to a directory.
* Open a """command window""", go to the directory with caseControlPower.py, type
-->`> caseControlPower.py`
->to start a graphical user interface, or use
-->`> caseControlPower.py --gui=False`
->to use interactive user input. You can use
-->`> caseControlPower.py --help`
->to check allowed command line options. You can also import this script as a module in another script and use the functions directly.

Please feel free to [ email me](mailto:bpeng@mdanderson.org ) if you have any question about this script.

### How to use this script
Our analysis assumes a disease predisposing locus and a marker locus that are
both in Hardy Weinberg Equilibrium. To calculate the power to detect a 
particular disease, you will need to disease mode (addition, recessive etc),
prevalence, disease allele frequency, marker allele frequency, linkage
disequilibrium between two loci (default to complete linkage), types of
cases and controls, genotype relative risk and significance level.

Attach:caseControlPower.png

Cases and controls should be specified as a dictionary of cases or controls with
different types. The pedigree types are  presented as a string of affection status
of father, mother, proband, and optional siblings of the proband. The affection
status is specified as 'A' for affected, 'U' for unaffected', and '*' for unknown.
For example, to calculate the power of two mixing pedigree types, you can use 
  cases = "{'**A': 100, '**AA': 200}"
which has 100 regular cases and 200 cases with an affected sibling. Of course, 
regular controls has the third letter as 'U'.

This script can carry out the following functions:

* '''Calculate power''':
->Given the number of cases of each pedigree type and genotype relative risk, calculate statistical power. You should specify number of cases and controls.

* '''Calculate minimal detectable relative risk'''
->Given number of cases and controls and power, calculate minimal detectable relative risk.

* '''Calculate needed samples from ratio between controls / cases'''
->Given power and relative risk, calculate number of cases and controls from a ratio.

* '''Calculate needed cased from fixed controls.'''
->Given power and relative risk, calculate number of cases and controls from a fixed number of controls.

* '''Calculate needed cased from fixed cases.'''
->Given power and relative risk, calculate number of cases and controls from a fixed number of cases.

This script accepts a large number of parameters and some of them are
hidden in the sense that they are not accessible from the graphical user
interface. Please refer the description of each parameter for details.

In addition to command line (you can use --gui=False to use this script in batch
mode) and a graphical user interface, you can also import this script and use 
the functions directly from another script. This will allow you to run a large
number of analyses easily and use a logging object to dump important internal
steps. Please refer to script [ PowerAnalysis.py](Attach:PowerAnalysis.py ) for an example. That script
loads caseControlPower.py and performs all calculates for Peng et al, 2010.

### Reference to command line parameters

python caseControlPower.py -h
A power calculator for case control association studies with know family
histories

This program is a power calculator for case-control association studies.
Compared to other power calculators, this program is unique in that it assumes
that we know the family history of cases and/or contorls so although only
samples are genotyped and used for statistical analysis, this study design has
higher power than regular case control studies because the probands have higher
probability to have disease alleles at the disease causing locus.  Cases and
controls should be specified as a dictionary of cases or controls with different
types. The pedigree types are presented as a string of affection status of
father, mother, proband, and optional siblings of the proband. The affection
status is specified as 'A' for affected, 'U' for unaffected', and '*' for
unknown. For example, to calculate the power of two mixing pedigree types, you
can use cases = "{'**A': 100, '**AA': 200}" which has 100 regular cases and 200
cases with an affected sibling. Of course, regular controls has the third letter
as 'U'.  This script can carry out the following functions:  1. Power
calculation: Given the number of cases of each pedigree type and genotype
relative risk, calculate statistical power. You should specify number of cases
and controls.  2. Calculate minimal detectable relative risk. Given number of
cases and controls and power, calculate minimal detectable relative risk.  3.
Calculate needed samples from ratio between ctrols / cases Given power and
relative risk, calculate number of cases and controls from a ratio.  4.
Calculate needed cased from fixed controls. Given power and relative risk,
calculate number of cases and controls from a fixed number of controls.  4.
Calculate needed cased from fixed cases. Given power and relative risk,
calculate number of cases and controls from a fixed number of cases.  This
program allows the users to variate the following passed-in parameters:  K:
Disease prevalence p: Disease allele frequency x: Marker allele frequency LD:
LD(D') alpha: Significant level power: Power Grr: Genotype relative risk  so
analyses on multiple disease models can be performed all at once.

usage: caseControlPower.py [--opt[=arg]] ...

options:
  -h, --help
        Display this help message and exit.

  --config=ARG (default: None)
        Load parameters from a configuration file ARG.

  --optimized
        Run the script using an optimized simuPOP module.

  --gui=[batch|interactive|True|Tkinter|wxPython] (default: None)
        Run the script in batch, interactive or GUI mode.

  --mode=ARG  (default: 'Additive')
        Disease models with different relationship between risks of AA, AB and
        BB, and different definitions of genotype relative risk. Let A be the
        disease allele, and r_AA, r_AB and r_BB be the penetrance of genotypes
        AA, AB and BB respectively, six disease types are provided for this
        analysis. Additive_AA and LogAdditive are alternatives of types Additive
        and Multiplicative by using alternative definitions of genotype relative
        risks and minimal detectable relative risk.
          Dominant:  r_AB=r_AA, Grr = r_AA / r_BB
          Recessive: r_AB=r_BB, Grr = r_AA / r_BB
          Multiplicative: r_AB=sqrt(r_AA * r_BB), Grr = r_AB / r_BB
          Additive:  r_AB=(r_AA + r_BB)/2, Grr = r_AB / r_BB
          Additive_AA: r_AB=(r_AA + r_BB)/2, Grr = r_AA / r_BB, Minimal
        Detectable Relative Risk = (p*p*Grr+pq(Grr+1))/(p*p+2pq)
          LogAdditive: r_AB=sqrt(r_AA * r_BB), Grr = r_AA / r_BB

  --K=ARG  (default: [0.05])
        Disease prevalence

  --p=ARG  (default: [0.15])
        Disease allele frequency

  --x=ARG  (default: [0])
        Marker allele frequency. If LD is 1, this parameter is ignored because
        it will be assumed to be the same as disease allele frequency.

  --genoP=ARG  (default: [0])
        This option is hidden. If specified, it overrides parameter p and
        specifies the frequency of causal genotypes so that
          p=sqrt(genoP) for a recessive model (frequency of AA).
          p=1-sqrt(1-genoP) for other models (frequency of AA + AB)

  --genoX=ARG  (default: [0])
        This option is hidden. If specified, it overrides paramter x and
        specifies the frequency of causal genotypes at the marker locus so that
          p=sqrt(genoX) for a recessive model (frequency of AA).
          p=1-sqrt(1-genoX) for other models (frequency of AA + AB)

  --Dprime=ARG  (default: 1)
        Linkage disequilibrium (D') between marker locus and disease locus. Only
        one of parameters Dprime and R2 should be specified.

  --R2=ARG  (default: 1)
        Linkage disequilibrium (R2) between marker locus and disease locus. Only
        one of the parameters Dprime and R2 should be specified.

  --analysis=ARG  (default: 'Statistical power')
        Perform different types of calculations, each require the specification
        of different parameters. More specifically,
        Statistical power: fix number of cases and controls, alpha, Grr,
        calculate power. Ignore ratio.
        Minimal detectable relative risk: fix number of cases and controls,
        alpha, power, calculate Grr. Ignore ratio.
        Sample size from #ctrl/#case: number of cases and controls treats as
        weights, fix ratio, alpha, Grr, power, calculate sample size.
        Sample size from fixed #ctrl: Fix number of controls, alpha, Grr, power
        and calculate number of cases. Ignore ratio. Numbers of different types
        of cases are treated as weights.
        Sample size from fixed #cases: Fix number of cases, alpha, Grr, power
        and calculate number of controls. Ignore ratio. Number of different
        types of controls are treated as weights.

  --cases=ARG  (default: {'**A': 1000})
        Type and number of cases. The count will be considered as weight in
        'sample size from #ctrl/#case', and 'sample size from fixed #ctrl'
        analysis

  --controls=ARG  (default: {'**U': 1000})
        Type and number of controls. The count will be considered as weight in
        'Sample size from #ctrl/#case'.

  --ratio=ARG  (default: 1)
        Ratio, use for calculate sample size from #ctrl/#case

  --alpha=ARG  (default: [1e-07])
        Significant level

  --power=ARG  (default: [0.8])
        Power of a case-control assocition test if number cases and relative
        risk are specified

  --Grr=ARG  (default: [1.2])
        Relative genotype risk for a disease model. The exact definition of this
        quantity depends on the disease model.
