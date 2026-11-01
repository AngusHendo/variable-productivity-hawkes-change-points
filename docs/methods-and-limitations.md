# Methods and limitations

## Computational methods

The source implements branching-process simulation for spatiotemporal Hawkes processes with exponential temporal and Gaussian spatial kernels. Productivity can be constant, piecewise constant or represented by several smooth functions.

The finite-window likelihood accounts for the temporal boundary and rectangular spatial boundary through Gaussian distribution-function terms. Shared parameters can be fitted under a stationary model before fitting segment-specific productivity values and candidate change times.

The stochastic-EM prototype treats the unknown change time as latent. It samples candidate times, weights them using the variable-productivity likelihood, updates the shared and segment-specific parameters, and narrows the candidate region over successive iterations.

## Simulation findings

The saved NeSI study contains 4,382 completed runs from a planned 4,500-task array. The main empirical pattern is that change-time recovery improves as the difference between pre- and post-change productivity increases. Similar segment productivities can leave the change time weakly identified even when the productivity estimates remain reasonably accurate.

Recovery also depends on the direction of change, the absolute productivity levels and the local allocation of events around the change. A larger total event count does not by itself guarantee precise change-time estimation.

## Limitations

These methods are research prototypes rather than a supported R package. The stochastic-EM standard-error calculation based on Louis's identity was not resolved before handover and is not included in the runnable example.

The Gambian disease analysis in the historical repository was preliminary. Case-level data and the unfinished multivariate extension are not included here. The later paper's Bayesian estimator, credible regions and vaccine-rollout interpretation are also outside the scope of this repository.
