#  Complete scripts for various applications

The following scripts were originally distributed with simuPOP in the `scripts` directory. They have been used for different research topics and can serve as examples on how to use simuPOP in real-world applications. Note that the documentations may be incorrect and there may be bugs in the script. %red%Please feel free to [contact me](mailto:bpeng@mdanderson.org) for any problem you have with these scripts.%% 


### Published papers that make use of simuPOP

These scripts were written for particular research papers. I have tried to update them so that they could run correctly under the current simuPOP version.

* [simuCDCV](Cookbook/simuCDCV) simulation for Reich's paper: On the allelic spectrum of human disease plus migration and complexity of the common disease. (for publication [ Peng and Kimmel, Genetics, 2007](simuCDCV](Cookbook/simuCDCV) simulation for Reich's paper: On the allelic spectrum of human disease plus migration and complexity of the common disease. (for publication [ Peng and Kimmel, Genetics, 2007)))

* [simuComplexDisease](Cookbook/simuComplexDisease) generate dataset for common complex disease  with certain number of disease susceptibility loci.  (for publication [ Peng et al, PLoS Genetics,  2007](simuComplexDisease](Cookbook/simuComplexDisease) generate dataset for common complex disease  with certain number of disease susceptibility loci.  (for publication [ Peng et al, PLoS Genetics,  2007)))

* [simuGWAS](Cookbook/simuGWAS) generate datasets by evolving HapMap samples forward in time. (for publication [Peng and Amos, BMC Bioinformatics, 2010](simuGWAS](Cookbook/simuGWAS) generate datasets by evolving HapMap samples forward in time. (for publication [Peng and Amos, BMC Bioinformatics, 2010)))

* [simuRareVariants](Cookbook/simuRareVariants) (SRV) simulate populations with rare variants. (described in [Peng and Liu, Hum Hered. 70(4): 287–291](simuRareVariants](Cookbook/simuRareVariants) (SRV) simulate populations with rare variants. (described in [Peng and Liu, Hum Hered. 70(4): 287–291))).

### Script for fun (research related)

These scripts are written in simuPOP for various purposes. They can be useful examples for your research.

* [landscapeGenetics](landscapeGenetics)) A landscape genetics simulation with different selection forces on X and Y axis and epistasis.

### Miscellaneous
* A power calculator for case control association studies with known family histories. Script [ caseControlPower.py]( Cookbook/caseControlPower ) calculates statistical power of case control association studies with known family histories. It makes use of simuPOP's simuOpt module for its graphical user interface and simuPOP.gsl module for probability functions. Please refer to Peng et al, Human Genetics, 2010 for details.
