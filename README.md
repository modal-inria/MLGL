# MLGL

[![Travis build status](https://travis-ci.com/modal-inria/MLGL.svg?branch=master)](https://travis-ci.com/modal-inria/MLGL) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/modal-inria/MLGL?branch=master&svg=true)](https://ci.appveyor.com/project/modal-inria/MLGL)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MLGL)](https://cran.r-project.org/package=MLGL) [![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/MLGL?color=blue)](http://cranlogs.r-pkg.org/badges/grand-total/MLGL) [![Downloads](https://cranlogs.r-pkg.org/badges/MLGL)](https://cran.rstudio.com/web/packages/MLGL/index.html)

The code was originally on an [R-forge repository](https://r-forge.r-project.org/projects/hcgglasso/).


This package implements a new procedure of variable selection in the context of redundancy between explanatory variables, which holds true with high dimensional data.


## Installation

From github:
``` r
library(devtools)
install_github("modal-inria/MLGL")
```

From CRAN:
``` r
install.packages("MLGL", repos = "https://cran.rstudio.com")
```

## Credits

**MLGL** is developed by Quentin Grimonprez, Guillemette Marot, Alain Celisse and Samuel Blanck.

Copyright Inria - Université de Lille

## Licence

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU Affero General Public License](https://www.gnu.org/licenses/agpl-3.0.en.html) for more details.


## References

*  Q. Grimonprez, S. Blanck, A. Celisse, G. Marot, MLGL: An R package implementing correlated variable selection by hierarchical clustering and group-Lasso. https://hal.inria.fr/hal-01857242

