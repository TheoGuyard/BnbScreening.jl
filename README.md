# Node-screening tests for the L0-penalized least-square problem

This repository contains numerical procedures linked to the paper

> Guyard, T., Herzet, C., & Elvira, C. (2022, May). Node-Screening Tests For The L0-Penalized Least-Squares Problem. In *ICASSP 2022-2022 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)* (pp. 5448-5452). IEEE.

available [here](https://arxiv.org/pdf/2110.07308.pdf). If you encounter a bug or something unexpected, please let me know by raising an issue on the project page or by contacting me by [mail](mailto:theo.guyard@insa-rennes.fr).

## Requirements

This repository works with Julia v1.5+. Please refer to the [Julia download page](https://julialang.org/downloads/) for install instructions.

**Dependencies :** [JuMP.jl](https://github.com/jump-dev/JuMP.jl), [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl), [Distributions.jl](https://github.com/JuliaStats/Distributions.jl), [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl)

Our code works with the [CPLEX](https://www.ibm.com/fr-fr/analytics/cplex-optimizer) solver through [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl). You cannot use CPLEX without having purchased and installed a copy of CPLEX Optimization Studio from [IBM](https://www.ibm.com). However, CPLEX is available for free to [academics and students](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students). We recommend you to follow carefully the install instructions described at [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl) to install the package.

## Install instructions

1. Clone the repository
```bash
git clone https://gitlab.insa-rennes.fr/Theo.Guyard/bnb-screening.git
```

2. Enter in the project root folder
```bash
cd bnb-screening
```

3. Ensure that Julia v1.5+ is installed
```bash
julia -v
> julia version 1.5.3 (or higher)
```

4. Install dependencies using [Julia's Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/) in REPL mode
```julia
(@v1.5) pkg> activate .
(bnb-screening) pkg> instantiate
```
It is possible that you also have to build [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl) manually with the command
```julia
(bnb-screening) pkg> build CPLEX
```

## ICASSP experiments

Experiments presented in the paper submitted to [ICASSP 2022](https://2022.ieeeicassp.org) are located in the folder `bnb-screening/exp/ICASSP/`. The available experiments are 
* `performances.jl` corresponding to results of table 1

Run experiments directly from the project root folder as follows.
```bash
julia --project=. -t 1 exp/ICASSP/performances.jl
```

The `--project=.` flag installs and builds dependencies if needed. The `-t 1` flag restricts computations to only 1 CPU core in order to avoid bias due to parallelization capabilities. Logs of the experiment are saved in the same folder as `performances.jl`. You can edit the experimental setup directly in the script, which is well commented.

## Running demo examples

The `bnb-screening/examples/` folder contains portions of code explaining how to use our package. You can run these scripts as follows.
```bash
julia --project=. examples/solve_bnb.jl
julia --project=. examples/solve_mip.jl
```

## Licence

This software is distributed under the [MIT Licence](https://mit-license.org).

## Cite this work

If you use this package for your own work, please consider citing it as :

```bibtex
@inproceedings{guyard2022node,
  title={Node-Screening Tests For The L0-Penalized Least-Squares Problem},
  author={Guyard, Th{\'e}o and Herzet, C{\'e}dric and Elvira, Cl{\'e}ment},
  booktitle={ICASSP 2022-2022 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)},
  pages={5448--5452},
  year={2022},
  organization={IEEE}
}
```
