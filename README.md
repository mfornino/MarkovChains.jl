# MarkovChains

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mfornino.github.io/MarkovChains.jl/)
<!--[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mfornino.github.io/MarkovChains.jl/dev)
[![Build Status](https://github.com/mfornino/MarkovChains.jl/workflows/CI/badge.svg)](https://github.com/mfornino/MarkovChains.jl/actions)
[![Coverage](https://codecov.io/gh/mfornino/MarkovChains.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mfornino/MarkovChains.jl)-->

## Introduction

A Julia package for Markov Chains and related functions. The roadmap for development includes:

* Defining a type hierarchy with Continuous and Discrete time MC further specialized to Univariate and Multivariate (in the spirit of [Distributions.jl](https://github.com/JuliaStats/Distributions.jl))
* Creating a number of convenience constructors. For instance, create a UnivariateContinuousTimeMarkovChain object from the discretization of a ItoDiffusionProcess object.
* Defining a method to compute the stationary distribution for each type
* Defining a method to extract a random sample from the MC
* Defining a method to determine whether the MC is irreducible

As of v0.1.3, I have developed the type ContinuousTimeMarkovChain, as well as the associated constructors and methods to compute the stationary distribution and drawing random samples.

## Resources

* **Documentation**: <https://mfornino.github.io/MarkovChains.jl/>

* **Support**: Please don't hesitate to contact me if you are interested in using the package.