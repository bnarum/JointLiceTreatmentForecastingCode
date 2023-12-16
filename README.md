# Joint Forecasting of Salmon Lice and Treatment Interventions in Aquaculture Operations

> by [Benjamin S. Narum](https://www.nhh.no/en/employees/faculty/benjamin-narum/) and [Geir D. Berentsen](https://www.nhh.no/en/employees/faculty/geir-drage-berentsen/)

This is code accompanying the paper "Joint Forecasting of Salmon Lice and Treatment Interventions in Aquaculture Operations". Estimation is done using the R-package [`TMB`](https://cran.r-project.org/web/packages/TMB/index.html). The negative log-likelihood functions are contained in the files `nll_lice.cpp` and  `nll_treat.cpp`, while the procedure to estimate the models can be found in the file `model_estimation.R`. The data is contained in the folder `data/`.
