
<!-- README.md is generated from README.Rmd. Please edit that file -->
![**funGp**](https://drive.google.com/thumbnail?id=1GsHPb5GnE0YS7b7nRpHqXyFKeHuLw3hF)

Gaussian process regression models with emphasis on functional inputs
---------------------------------------------------------------------

**Description:** funGp is a regresssion package based on Gaussian process models. It allows inputs to be either scalar, functional (represented as vectors), or a combination of both. A dimension reduction functionality is implemented in order to reduce the complexity of the functional inputs while keeping a sufficient amount of information to use in the model. Moreover, funGp offers a model selection feature which allows to tune different characteristics of the model such as the active scalar and functional inputs, the type of kernel function and the family of basis function used for the projection of the inputs. This is an extension of the work presented in Betancourt et al. (2020).

**Main functionalities** <br /> :small\_blue\_diamond: Creation of regression models <br /> :small\_blue\_diamond: Output estimation at unobserved input points <br /> :small\_blue\_diamond: Random sampling from a Gaussian process model <br /> :small\_blue\_diamond: Heuristic optimization of model structure <br />

**Note:** funG pwas first developed in the frame of the RISCOPE research project, funded by the French Agence Nationale de la Recherche (ANR) for the period 2017-2021 (ANR project No. 16CE04-0011, [RISCOPE.fr](https://perso.math.univ-toulouse.fr/riscope/)), and certified by SAFE Cluster.

This project is licensed under the GPL-3 License. <br /><br />

**Installation**

    # Install release version from CRAN
    install.packages("funGp")

    # Install development version from GitHub
    # way 1
    library(devtools)
    install_github("djbetancourt-gh/funGp", dependencies = TRUE)

    # way 2
    library(githubinstall)
    githubinstall("funGp", dependencies = TRUE)

<br />

**Manual** :book: <br /> [Gaussian Process Regression for Scalar andFunctional Inputs with funGp - The in-depth tour](https://drive.google.com/file/d/1MtYi-Qq-BNZpbp1SWWG4Fb35HvVqRHQM/view?usp=sharing) <br /><br />

**Authors:** José Betancourt :wrench: (IMT, ENAC), François Bachoc (IMT) and Thierry Klein (IMT, ENAC).

**Contributors:** Déborah Idier (BRGM) and Jérémy Rohmer (BRGM).

:wrench: maintainer - [jbetanco@math.univ-toulouse.fr](jbetanco@math.univ-toulouse.fr) <br /><br />

**Acknowledgments:** we are grateful to Yves Deville from Alpestat for his advice on the documentation of R packages. <br /><br />

**References** <br />

Betancourt, J., Bachoc, F., Klein, T., Idier, D., Pedreros, R., and Rohmer, J. (2020), "Gaussian process metamodeling of functional-input code for coastal flood hazard assessment". *Reliability Engineering & System Safety*, **198**, 106870. [RESS](https://www.sciencedirect.com/science/article/abs/pii/S0951832019301693) - [HAL](https://hal.archives-ouvertes.fr/hal-01998724)

Betancourt, J., Bachoc, F., Klein, T., and Gamboa, F. (2020), Technical Report: "Ant Colony Based Model Selection for Functional-Input Gaussian Process Regression. Ref. B3D-WP3.2". *RISCOPE project*. [HAL](https://drive.google.com/file/d/1GnalLS9jEr9AxPKmQk0S1bLQ7whuLm1T/view?usp=sharing)

Betancourt, J., Bachoc, F., and Klein, T. (2020), R Package Manual: "Gaussian Process Regression for Scalar and Functional Inputs with funGp - The in-depth tour". *RISCOPE project*. [HAL](https://drive.google.com/file/d/1MtYi-Qq-BNZpbp1SWWG4Fb35HvVqRHQM/view?usp=sharing)
