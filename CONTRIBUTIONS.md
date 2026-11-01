# Contribution and provenance

Angus Henderson developed the early R implementation collected here while working with Dr Conor Kresin at the University of Otago from September 2024 to February 2025. The formal Assistant Research Fellow contract ran from November 2024 to February 2025.

The work included:

- univariate and spatiotemporal Hawkes-process simulation;
- finite-window likelihood and compensator calculations;
- maximum-likelihood fitting for constant and time-varying productivity;
- single- and multiple-change-point prototypes;
- a stochastic-EM prototype treating the change time as latent;
- profile-likelihood and bootstrap uncertainty investigations; and
- large simulation studies on New Zealand's NeSI infrastructure.

Conor Kresin provided research direction and helped derive and correct parts of the likelihood. The implementation should not be read as a claim that every underlying mathematical result originated independently with Angus.

Sophie Phillips subsequently continued the research, including separate SEM and two-stage developments in the historical repository. Those later folders are not included here. The resulting paper by Conor Kresin, Boris Baeumer and Sophie Phillips acknowledges Angus for developing a stochastic EM algorithm for change-point detection.

This repository is a curated computational record. It is not the implementation of the paper's later Bayesian MAP and centre-of-mass method.
