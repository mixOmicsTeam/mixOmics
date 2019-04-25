# mixOmics


This repository includes the package now hosted on BioConductor and our current development version!


**Installation**
* Step 0. Mac OS users only: install X Quartz first https://www.xquartz.org/

* Step 1. To obtain the latest update of mixOmics, you will need to pull from our gitHub page via the devtools and the install_github libraries. Install the libraries 'devtools' in R, then load and install the latest stable version of mixOmics from gitHub (as bug-free as it can be):

> library('devtools')

> install_github("mixOmicsTeam/mixOmics")

You can also install the development version  
> install_github("mixOmicsTeam/mixOmics", ref="devel")

Check after install that the following does not throw any error (see step 0) and that the welcome message confirms you have installed version **>** 6.7.1:

>library(mixOmics) 
Loaded mixOmics 6.7.1


**Bug reports**
If you would like to report a bug or issue: https://github.com/mixOmicsTeam/mixOmics/issues

Thank you for using mixOmics!

**About the team**
mixOmics is collaborative project developed by the mixOmics team (Kim-Anh Lê Cao - https://lecao-lab.science.unimelb.edu.au, Florian Rohart - http://florian.rohart.free.fr, Ignacio González and Sébastien Déjean - https://perso.math.univ-toulouse.fr/dejean/), key contributors (Benoît Gautier, François Bartolo) and several key collaborators. The project started at the Institut de Mathématiques de Toulouse, Université Paul Sabatier, Toulouse, France and was then further extended in Australia, at the University of Queensland, Brisbane (2009 – 2016) and at the University of Melbourne, Australia (2017 – ). We focus on statistical methods development for biological data integration and implementation in R.


**Why this toolkit?**
mixOmics offers a wide range of novel multivariate methods for the exploration and integration of biological datasets with a particular focus on variable selection. Single ‘omics analysis does not provide enough information to give a deep understanding of a biological system, but we can obtain a more holistic view of a system by combining multiple ‘omics analyses. Our mixOmics R package proposes a whole range of multivariate methods that we developed and validated on many biological studies to gain more insight into ‘omics biological studies.


**Want to know more?**

www.mixOmics.org (tutorials and resources)

https://mixomicsteam.github.io/Bookdown/ Our latest bookdown vignette.



