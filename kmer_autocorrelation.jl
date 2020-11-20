### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 42cb734c-2b47-11eb-12de-6d2181461d88
begin
	using LinearAlgebra
	using BioSequences
	using Plots
	using Random
	using Statistics
end

# ╔═╡ 52fdcf9e-2b47-11eb-3710-6f1ef98c4e51
const KERNEL = let
	equations = []
	
	# Constrain one: All scaled frequencies sum to zero
	push!(equations, ones(Int, 256))
	
	# Constrain two: sum(xABC) == sum(ABCx)
	for base in 0:63
		equation = zeros(Int, 256)
		for affix in 0:3
			equation[base << 2 + affix + 1] += 1
			equation[affix << 6 + base + 1] -= 1
		end
		push!(equations, equation)
	end
	
	# Return the nullspace
	nullspace(Matrix(transpose(reduce(hcat, equations))))
end;

# ╔═╡ fea5f9c0-2b47-11eb-268b-03414fc4887e
"Use the kernel to project kmer freq matrix to orthonormal matrix"
function project_kmer_freq(freq::AbstractMatrix{<:Number})
	mean = sum(freq, dims=2) / size(freq, 2)
	shifted = freq .- mean
	return shifted * KERNEL
end;

# ╔═╡ 7e06120e-2b48-11eb-10fa-5765d37c680a
"Generate kmer spectrum from random nucleotides"
function random_kmer_spectrum(len::Integer, N::Integer)
	seq = LongSequence{DNAAlphabet{2}}(len)
	factor = 1 / (length(seq) - 3)
	counts = zeros(Int, 256)
	result = zeros(Float64, (N, 256))
	
	for row in 1:N
		rand!(seq.data)
		fill!(counts, zero(eltype(counts)))
	
		@inbounds for mer in each(DNAMer{4}, seq)
			index = reinterpret(Int, mer.fw) + 1
			counts[index] += 1
		end

		@inbounds for i in eachindex(counts)
			result[row, i] = counts[i] * factor
		end
	end
	
	return result
end

# ╔═╡ 9609e20c-2b4a-11eb-0f01-016bc09728d0
"Obs x features matrix"
function auto_correlation(matrix::AbstractMatrix)
	correlation = cor(matrix, dims=1)
	
	# Zero out correlation with a feature and itself, so as not
	# to saturate the heatplots below needlessly
	for i in 1:size(correlation, 1)
		correlation[i, i] = 0
	end
	correlation
end

# ╔═╡ 605618d0-2b4d-11eb-1bd3-ad2a83c52183
begin
	spectrum = random_kmer_spectrum(2_000, 10_000)
	orthonormal = project_kmer_freq(spectrum)
	corr_spectrum = auto_correlation(spectrum)
	corr_orthonormal = auto_correlation(orthonormal)
	
	mn = min(minimum(corr_spectrum), minimum(corr_orthonormal))
	mx = max(maximum(corr_spectrum), maximum(corr_orthonormal))
	clim = (mn, mx)
end;

# ╔═╡ d424aef0-2b4a-11eb-1e20-d93067c6038a
heatmap(corr_spectrum, clim=clim)

# ╔═╡ 70b600e0-2b49-11eb-1229-5fc87c84fc1e
heatmap(corr_orthonormal, clim=clim)

# ╔═╡ Cell order:
# ╠═42cb734c-2b47-11eb-12de-6d2181461d88
# ╠═52fdcf9e-2b47-11eb-3710-6f1ef98c4e51
# ╠═fea5f9c0-2b47-11eb-268b-03414fc4887e
# ╠═7e06120e-2b48-11eb-10fa-5765d37c680a
# ╠═9609e20c-2b4a-11eb-0f01-016bc09728d0
# ╠═605618d0-2b4d-11eb-1bd3-ad2a83c52183
# ╠═d424aef0-2b4a-11eb-1e20-d93067c6038a
# ╠═70b600e0-2b49-11eb-1229-5fc87c84fc1e
