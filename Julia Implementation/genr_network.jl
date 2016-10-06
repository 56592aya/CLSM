__precompile__
module GenrNetwork
export N ,K ,adj_matrix, η0, η1, α, ϵ, Θ, β
using Distributions
using GenericFuncs:sort_by_argmax

N=100
K=13
ϵ = 1e-30
η0 = 10
η1 = 1

β_diag = zeros(Float64, (K,K))

dbeta=Beta(η0, η1)
β_kk = rand(dbeta, K)
for k in 1:K
    β_diag[k,k] = β_kk[k]
end
β=β_diag + ϵ*(ones(Int64, (K,K))-eye(Int64, K))


α=Array(Float64, K)
for k in 1:K
    α[k] = 0.01
end

ddir = Dirichlet(α)
Θ = zeros(Float64, (N,K))
for a in 1:N
    Θ[a,:] = rand(ddir)
end
Θ=sort_by_argmax(Θ)

adj_matrix = zeros(Int64, (N,N))
for a in 1:N
    for b in 1:N
        if(b > a)
            dmult_a = Multinomial(1,vec(Θ[a,:]))
            dmult_b = Multinomial(1,vec(Θ[b,:]))
            z_send_ab = rand(dmult_a,1)
            z_recv_ab = rand(dmult_b,1)
            z_send_ab_idx = indmax(z_send_ab)
            z_recv_ab_idx = indmax(z_recv_ab)
            dbinom = Binomial(1,β[z_send_ab_idx, z_recv_ab_idx])
            adj_matrix[a,b]=adj_matrix[b,a] = rand(dbinom, 1)[1]
        end
    end
end
end
