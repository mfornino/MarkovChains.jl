"""
ContinuousTimeMarkovChain{T}

Type for storing a continuous time Markov chain.

Can be constructed directly by supplying the infinitesimal generator infinitesimal_generator, or by supplying an ItoDiffusionProcess to be approximated.

See also https://en.wikipedia.org/wiki/Continuous-time_Markov_chain#Stationary_distribution

"""
struct ContinuousTimeMarkovChain{T}
	
	infinitesimal_generator::AbstractMatrix{<:Real}
	states::AbstractVector{T}

	# Standard Constructor with Consistency Checks
	function ContinuousTimeMarkovChain(infinitesimal_generator::AbstractMatrix{<:Real}, 
									   states::AbstractVector{T}) where T 
		
		# Check infinitesimal generator is a square matrix
		if ~ size(infinitesimal_generator, 1) == size(infinitesimal_generator, 2)
			error("infinitesimal generator must be a squadre matrix.")
		end
		
		# Check infinitesimal generator rows sum to 0
		# Numerically, this means that the distance between sum(row_i) and 0.0 
		# can at most be the machine precision around the value of the diagonal 
		# element. On top of that we add 2 orders of magnitude.
		row_sums = abs.(sum(infinitesimal_generator, dims=2))
		T_i = typeof(infinitesimal_generator[1,1])
		if T_i <: Integer
			if any(row_sums .> zero(T_i))
				error("Rows of the infinitesimal generator must sum to 0.")
			end
		else
			check = false
			for (idx, row_sum) in enumerate(row_sums)
				check = ~isapprox(row_sum, zero(T_i), atol=1e2 .* eps(infinitesimal_generator[idx,idx]))
				check ? error("Rows of the infinitesimal generator must sum to 0.") : continue
			end
		end

		# Check sign restriction
		if any(diag(infinitesimal_generator, 0) .> 0)
			error("Diagonal elements of the infinitesimal generator must be negative or zero.")
		end

		# Check infinitesimal generator and supplied state vector have compatible sizes
		if size(infinitesimal_generator, 1) != length(states)
			error("Size mismatch between infinitesimal generator infinitesimal_generator and the state vector states.")
		end

		# Allocate object
		new{T}(infinitesimal_generator, states)
	end
	
end


# Constructor which omits State Vector, assume state space is a range vector.
function ContinuousTimeMarkovChain(infinitesimal_generator::AbstractMatrix{<:Real}) 

	# If State Vector not provided, then use 1:size(infinitesimal_generator,1) as placeholder
	states = Vector{Int64}(1:size(infinitesimal_generator, 1))

	# Allocate object
	return ContinuousTimeMarkovChain(infinitesimal_generator, states)
end


# Constructor for single-dimensional Diffusion Process	
function ContinuousTimeMarkovChain(dp::ItoDiffusionProcess, grid::Vector{<:Real})

	# Grab length of State Vector
	n = length(grid)
	
	# Preallocate intensity vectors
	diag = Vector{Float64}(undef, n)
	fwd = similar(diag[1:end-1])
	bwd = similar(fwd)
	dS = Vector{Float64}(undef, n-1)

	# Unequally spaced grids are supported
	dS = diff(grid)
	dS_dwn = vcat(dS[1], dS)
	dS_up = vcat(dS, dS[end])
	dS_avg = (dS_up + dS_dwn) ./ 2
	# dS_avg[1] = dS_avg[1] ./2
	# dS_avg[end] = dS_avg[end] ./2

	# Compute drift and diffusion terms
	drift = dp.μ.(grid)
	diffusion = dp.σ.(grid) .^ 2 ./ 2

	# Compute intensity up and down using Upwind Scheme
	bwd = - min.(drift, 0) ./ dS_dwn + diffusion ./ (dS_avg .* dS_dwn)
	fwd = max.(drift, 0) ./ dS_up + diffusion ./ (dS_avg .* dS_up)

	# By construction, outflow intensity on diagonal s.t. rows sum to 0
	bwd[1], fwd[end] = 0, 0
	diag = - fwd - bwd

	# Populate Infinitesimal Generator as a Sparse Matrix
	infinitesimal_generator = spdiagm(-1 => bwd[2:end], 0 => diag, 1 => fwd[1:end-1])

	# Allocate Object
	return ContinuousTimeMarkovChain(infinitesimal_generator, grid)
end


# Overloading the show operator to get pretty-printing behavior.
function Base.show(io::IO, ::MIME"text/plain", mc::ContinuousTimeMarkovChain) 
	print(io, "Continuous Time Markov Chain object.\n")
	bigchain = length(mc.states) > 10
	if bigchain
		print(io, "Stored objects are too large. Showing only first few elements.\n")
		print(io, "\nInfinitesimal generator:\n")
		display(Matrix(mc.infinitesimal_generator[1:10, 1:10]))
		print(io, "\nState space:\n")
		display(mc.states[1:10])
	else
		print(io, "\nInfinitesimal generator:\n")
		display(Matrix(mc.infinitesimal_generator))
		print(io, "\nState space:\n")
		display(mc.states)
	end
end


"""
	stationary_distribution(mc::ContinuousTimeMarkovChain, tol::T)

Computes the stationary distribution of the continuous time Markov chain mc.

By convention, the stationary distribution is treated as a PMF.

See also https://en.wikipedia.org/wiki/Continuous-time_Markov_chain#Stationary_distribution
"""
function stationary_distribution(mc::ContinuousTimeMarkovChain{T},
								 step_size::Real=1e8,
								 maxit::Integer=20,
								 tol::Real=1e-8) where T
	
	# Delta parameter for implicit method
	Δ = step_size

	# Initial guess is uniform over all states
	g0 = ones(Float64, size(mc.states)) ./ length(mc.states)
	g = similar(g0)
	
	# Compute RHS matrix in KFE
	B = I - Δ .* transpose(mc.infinitesimal_generator)

	# 
	it = 0
	supnorm = 1.0e6
	while supnorm > tol && it < maxit
		g = B \ g0
		supnorm = maximum(abs.(g .- g0)) 
		g0 = g
		it += 1
	end

	if it == maxit
		error("Algorithm for finding the stationary distribution did not converge.")
	end

	g ./= sum(g)

	if isa(Real, T)
		return DiscreteNonParametric(mc.states, g)
	else
		return Categorical(g)
	end
end




"""
	random_sample(mc::ContinuousTimeMarkovChain, Is::Integer=1, draws::Real=1000)

Given a starting index Is of type Integer, it returns a sequence of draws of length N.

The output Tuple contains a cumulative elapsed time T and the pseudorandom sequence of indices.

"""
function random_sample(mc::ContinuousTimeMarkovChain, Is::Integer=1, N::Real=1000)

	spy = pattern(mc)
	degenerate = size(spy,1) == 1

	traj = Vector{Int64}(undef, N + 1)
	traj[1] = Is

	T = Vector{Float64}(undef, N + 1)
	T[1] = 0

	if ~ degenerate
		for draw = 1:N

			Δ, traj[draw+1] = conditional_draw(mc, traj[draw], spy)

			T[draw+1] = T[draw] + Δ

		end
	else

		traj[2:end] = Is
		T[2:end] = Inf
	
	end

	return T, traj
end

function conditional_draw(mc::ContinuousTimeMarkovChain, orig_idx::Integer, pattern::AbstractMatrix{Tv}) where Tv<:Bool

	ran = 1:length(mc.states)

	cols = ran[view(pattern, orig_idx, :)]
	
	Δ, idx = findmin(rand.(Exponential.(1 ./ view(mc.infinitesimal_generator, orig_idx, cols))))

	return Δ, cols[idx]

end

# Better not to have to use this one!!!
function conditional_draw(mc::ContinuousTimeMarkovChain, orig_idx::Integer)

	return conditional_draw(mc, orig_idx, pattern(mc))

end

function pattern(mc::ContinuousTimeMarkovChain)

	return mc.infinitesimal_generator .> 0.0

end

