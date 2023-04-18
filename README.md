[PGSToolKit](https://github.com/swvanderlaan/PGSToolKit)<img align="right" height="250" src=LOGOS/fulllogo.png>
============
[![DOI](https://zenodo.org/badge/128260489.svg)](https://zenodo.org/badge/latestdoi/128260489) 

[![Languages](https://skillicons.dev/icons?i=bash,r,py)](https://skillicons.dev) 

A Toolkit to calculate Polygenic Scores (PGS) using _PLINK2_, _PRSice2_, _RapidoGS_, or _PRS-CS_. 

The **PGSToolKit** acts as a wrapper tool for several PGS methods and can be used to compute polygenic scores for a target population using genome-wide association study (GWSA) summary statistics and a reference (_e.g._ 1000G phase 3). Current methods for computing polygenic scores require a wide range of different input formats and parameters. **PGSToolKit** provides a streamlined and clean alternative by using a fixed data format and parameters, configurable from a single configuration file. In short, a configuration file, `pgstoolkit.conf`, is set, and `pgstoolkit.sh` can be readily submitted to the server. 

Current supported features are:

- Perform _quality control_ by filtering variants based on imputation score and minor allele frequency.
- Compute polygenic scores using _[PRS-CS](https://github.com/getian107/PRScs)_, _[RapidoPGS](https://github.com/GRealesM/RapidoPGS)_ and _[PRSice2](https://choishingwan.github.io/PRSice/)_ based on GWAS summary statistics.
- Compute polygenic scores using the [allelic score](https://www.cog-genomics.org/plink/2.0/score) function of _[PLINK2](https://www.cog-genomics.org/plink/2.0/)_ based on posterior variant effect sizes or weights.

#### Instructions and usage
Detailed instructions on using **PGSToolKit** can be found in the [Wiki](https://github.com/swvanderlaan/PRSToolKit/wiki).

#### Requirements
All scripts are annotated for debugging purposes - and future reference. **PGSToolKit** was tested within the context of a CentOS7 system with SLURM. 


--------------

#### The MIT License (MIT)
[Copyright (c)](copyright.md) 2017-2023 Sander W. van der Laan | s.w.vanderlaan [at] gmail [dot] com | [vanderlaanand.science](https://vanderlaanand.science) & [Anton Ligterink](https://www.linkedin.com/in/anton-ligterink/).
