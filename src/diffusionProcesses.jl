"""
ItoDiffusionProcess

Type for storing an Itô diffusion process. 

It consists of a function μ(x) for the drift term, and a function σ(x) for the diffusion term.

See also https://en.wikipedia.org/wiki/It%C3%B4_diffusion
"""
struct ItoDiffusionProcess
	μ::Function
	σ::Function
	# interval::Tuple{T, T} where T<:Real
	# reflective::Tuple{Bool, Bool}

	# function ItoDiffusionProcess(mu::Function, sigma::Function, interval::Tuple{T, T} where T<:Real)

	# 	# If reflective nature of barriers is not specified, assume true
	# 	reflective = (True, True)

	# 	# Allocate object
	# 	new(mu, sigma, interval, reflective)
	# end

end