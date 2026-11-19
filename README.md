# QuasiSeq

Analyzes RNA sequencing count tables with quasi-likelihood methods.

Compatibility updates by Renova.

Uses the supported R_Calloc and R_Free interfaces, updates two help pages, and removes the stale Date field.

## Installation

```r
install.packages("QuasiSeq", repos = c(
  Renova = "https://renova-r.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))
```

Version `1.0.11.0.9001` is based on [CRAN 1.0-11-0](https://cran.r-project.org/package=QuasiSeq).

## Credit

Steve Lund [aut, cre], Long Qu [aut, ctb], Klirk Bhasavanich [aut], Ian Marschner [aut] (The author of glm2::glm.fit2, which was modified slightly leading to glm.fit3 in this package.), R Core Team [aut] (The author of stats::glm.fit, which was modified slightly leading to glm.fit3 in this package.).

Licensed under `GPL (>= 2)`.
