# dlmt: Dynamic Linear Models with Transformations

Bare-bones package for fitting dynamic linear model with optimized outcome transformations, 
based on "Athlete rating in multi-competitor games with scored outcomes via monotone transformations" 
(Che and Glickman 2023+, [arxiv link](https://arxiv.org/abs/2205.10746)).


## Installation

The `dlmt` package uses `rstan`, which may require a bit of setup.
If you do not have the `rstan` package installed, see [this link](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) to get started.

Then `dlmt` can be installed using:
```
devtools::install_github("jche/dlmt")
```

## Usage

`dlmt` provides two functions:

1. `sim_data()`: simulate a sample dataset
2. `dlmt()`: fit dynamic linear model with transformations to a dataset

A simple use-case would be:
```
df <- sim_data()
results <- dlmt(df)
```

See `?sim_data` and `?dlmt` for further details about these functions and their outputs.

