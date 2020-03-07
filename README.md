# 2020-ncov

Analysis of the 2019-nCoV outbreak during 2019/20. _Note: this is working repository, so code and data are likely to change over time_

### Guide to files for `stoch_model`

This is a stochastic SEIR model implemented using Euler-Maruyama, with likelihood estimated using SMC by jointly fitting to cases in Wuhan and exported cases over time in countries with high connectivity to Wuhan.

Data loading and model run script is in `scripts/main_model.r`. Calls the following R files:

> `R/load_timeseries.r` - Load and format timeseries

> `R/model_functions.r` - Load process model and SMC

> `R/plotting_functions.r` - Plotting functions

> `outputs_main.R` - Run main model outputs

The code and data used for V1 of [our pre-print on early transmission dynamics](https://www.medrxiv.org/content/10.1101/2020.01.31.20019901v1) can be found in `stoch_model_V1_paper`, with same paths as above.

The code and data used for our final Lancet Infectious Diseases paper can be found in `stoch_model_V2_paper`, with same paths as above.


#### Reference

[Kucharski AJ, Russell TW, Diamond C et al. Early dynamics of transmission and control of 2019-nCoV: a mathematical modelling study. Lancet Infectious Diseases, 2020](https://www.medrxiv.org/content/10.1101/2020.01.31.20019901v1)