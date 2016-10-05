module GlobalInit
using GenrNetwork:K, N, adj_matrix, η0, η1
using Distributions

M = convert(Int64, sum(adj_matrix/2))
γ = zeros(Float64, (N, K))
for a in 1:N
    for k in 1:K
        dunif = Uniform(0,1)
        γ[a,k] = rand(dunif)
    end
end

total_pairs = convert(Int64, N*(N-1)/2)
ones_prob = M/total_pairs
τ0 = zeros(Float64, K)
τ1 = zeros(Float64, K)
for k in 1:K
    τ0[k] = η0
    τ1[k] = η1
end
export γ, τ0, τ1

end
