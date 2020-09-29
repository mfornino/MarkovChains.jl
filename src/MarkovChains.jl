module MarkovChains

include("diffusionProcesses.jl")
include("CTMC.jl")

using LinearAlgebra
using SparseArrays
using Distributions

import Base: 
	show

# Export types
export 
	ContinuousTimeMarkovChain,
	ItoDiffusionProcess

# Export methods
export 
	stationary_distribution,
	random_sample

"""
A Julia package for simulating Markov Chains.

"""


end
