<p align="center">
  <img src="https://raw.githubusercontent.com/tlcaputi/dlm/main/docs/assets/logo.svg" alt="dlm logo" width="120">
</p>

<h1 align="center">dlm</h1>

<p align="center">
  <strong>Distributed Lag Models for Stata</strong><br>
  Generalize event studies to continuous treatments
</p>

<p align="center">
  <a href="https://tlcaputi.github.io/dlm/">Documentation</a> &middot;
  <a href="https://github.com/tlcaputi/dlm">R version</a>
</p>

---

Stata implementation of the distributed lag model (DLM) framework from [Schmidheiny and Siegloch (2023, *Journal of Applied Econometrics*)](https://doi.org/10.1002/jae.2971).

DLMs generalize the canonical event study to settings with **continuous treatments** that can change in magnitude, sign, and timing throughout the study period. When applied to a binary absorbing treatment, the DLM produces estimates that are numerically identical to a binned-endpoint event study.

## Installation

```stata
net install dlm, from("https://raw.githubusercontent.com/tlcaputi/dlm-stata/main/")
```

Requires [reghdfe](http://scorreia.com/software/reghdfe/): `ssc install reghdfe`

## Quick Start

```stata
* Generate test data (true treatment effect = -3)
dlm_gen_data, n_groups(500) n_times(20) treat_prob(0.4) seed(42) clear

* Estimate DLM
dlm outcome, exposure(post) unit(unit) time(time) from(-3) to(3)

* View results
matrix list e(betas)
```

## Why DLMs?

The canonical event study uses event-time dummies, which only work for binary treatments that turn on once and stay on. Many empirical settings involve **continuous treatments** — tax rates, minimum wages, policy dosages — where event-time dummies don't apply. The DLM replaces these dummies with leads and lags of the treatment variable itself, naturally handling treatments that are continuous, change sign, vary in magnitude, and occur multiple times per unit.

See the [documentation](https://tlcaputi.github.io/dlm/) for a full explanation with concrete examples.

## Citation

> Schmidheiny, K. and S. Siegloch (2023). "On event studies and distributed-lag models: Equivalence, generalization and practical implications." *Journal of Applied Econometrics*, 38(5): 695-713.
