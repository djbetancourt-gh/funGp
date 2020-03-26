# funGp

## Gaussian process regression models with emphasis on functional inputs

**Description:** funGp is a regresssion package based on Gaussian process models. It allows inputs to be either scalar, functional (represented as vectors), or a combination of both. A dimension reduction functionality is implemented in order to reduce the complexity of the functional inputs while keeping a sufficient amount of information to use in the model. Moreover, funGp offers a model selection feature which allows to tune different characteristics of the model such as the active scalar and functional inputs, the type of kernel function and the family of basis function used for the projection of the inputs.

**Main functionalities** <br />
:small_blue_diamond: Creation of regression models <br />
:small_blue_diamond: Output estimation at unobserved input points <br />
:small_blue_diamond: Random sampling from a Gaussian process model <br />
:small_blue_diamond: Heuristic optimization of model structure <br />

**Note:** funG pwas first developed in the frame of the RISCOPE research project, funded by the French Agence Nationale de la Recherche (ANR) for the period 2017-2021 (ANR project No. 16CE04-0011, [RISCOPE.fr](https://perso.math.univ-toulouse.fr/riscope/)), and certified by SAFE Cluster.

This project is licensed under the GPL-3 License.
<br /><br />

**Installation**

```
# Install release version from CRAN
install.packages("funGp")
```

```
# Install development version from GitHub
# way 1
library(devtools)
install_github("djbetancourt-gh/funGp")

# way 2
library(githubinstall)
githubinstall("funGp")
```
<br />

**Authors:** José Betancourt :wrench: (IMT, ENAC), François Bachoc (IMT) and Thierry Klein (IMT, ENAC).

**Contributors:** Déborah Idier (BRGM) and Jérémy Rohmer (BRGM).

:wrench: maintainer - [jbetanco@math.univ-toulouse.fr](jbetanco@math.univ-toulouse.fr)
<br /><br />

**Acknowledgments:** we are grateful to Yves Deville from Alpestat for his advice on the documentation of R packages.
