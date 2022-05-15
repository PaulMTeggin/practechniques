# practechniques - Practical Actuarial Techniques

A companion R package for the online book Practical Actuarial Techniques (in development).

**THIS IS A PRE-RELEASE VERSION AND NO RELIANCE WHATSOEVER MAY BE PLACED ON IT.**

At the time of writing (May 2022) the main feature is a framework
using Monte Carlo simulation to calculate p-values for statistical goodness-of-fit tests
where the parameters have been estimated from the data.

This is a generalisation of the approach in the
[`KScorrect`](https://cran.r-project.org/web/packages/KScorrect/index.html) package
and has been tested against results from `KScorrect`.

The approach has been generalised in three ways:

* General distributions are supported, rather than the closed list used by `KScorrect`.
* Multiple statistical tests are supported, not just [KS](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test).
* Testing can be performed against distributions fitted to overlapping data, not just IID data,
using the idea of a Gaussian copula to induce autocorrelation consistent with overlapping data suggested in section 4.2 of the
[2019 paper](https://www.cambridge.org/core/journals/british-actuarial-journal/article/calibration-of-var-models-with-overlapping-data/B20D66D81DB918AFD3BBDF9EDAC20863)
by the [Extreme Events Working Party](https://www.actuaries.org.uk/practice-areas/life/research-working-parties/extreme-events)
of the UK [Institute and Faculty of Actuaries](https://www.actuaries.org.uk/).

With careful configuration, the framework here can in principle also be used where parameters are known in advance
rather than estimated from the data, but there is very little value to this use case, as Monte Carlo is rarely necessary
when the parameters are known (and is certainly not necessary for the KS and AD tests).

## Installation

There is as yet no official release.
The code in this repository is the development version,
which can be installed from github using the `remotes` package:

```R
require(remotes)
remotes::install_github('PaulMTeggin/practechniques')
```

If you don't have `remotes` installed you should first run

```R
install.packages('remotes')
```

## Bug reports 

Users of `practechniques` are encouraged to report bugs here 
(go to *issues* in the menu above, 
and press *new issue* to start a new bug report
or feature request).
