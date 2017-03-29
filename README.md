# Comparison different MCMC methods

Cool Walking and methods like J-Walking, S-Walking, Smart Darting, and Parallel Tempering are compared in [1]. In this repository we do not a full comparison (yet), but compare it with a new MCMC approach in which we do not used chains on different temperatures, but pick from a set of (symmetric) proposal distributions.

## Potential function

The potential function at different temperatures. 

![Potential function](https://github.com/mrquincle/cool-walking/blob/master/pictures/potential_energy.png?raw=true "Potential function")

For `T=0.1` we have only a few peaks. For `T=3.0` we have a much more uniform distribution.

The warmer distribution is used to initiate larger jumps.

## Results

The current results with adding complexity in the proposal distribution rather than smoothing the target distribution:

![MCMC hybrid](https://github.com/mrquincle/cool-walking/blob/master/pictures/mcmc_hybrids.png?raw=true "MCMC hybrid")

# References

[1] Cool Walking: A New Markov Chain Monte Carlo Sampling Method (Brown, Head-Gordon, 2002)
