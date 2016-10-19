__precompile__
module Inference
export ELBO

using DataFrames
using Distributions
using Gadfly
using GenericFuncs:get_links, get_nonlinks, get_static_neighbors, get_static_nonneighbors,Elog_Dirichlet,Elog_Beta, log_sum
using GenrNetwork



M = convert(Int64, sum(adj_matrix/2))
links=get_links(adj_matrix)
immutable MyNodePair
    source::UInt64
    sink::UInt64
end

non_links = get_nonlinks(adj_matrix)
neighbors = get_static_neighbors(adj_matrix)
non_neighbors = get_static_nonneighbors(adj_matrix)
deg_a = zeros(Int64, N)
for a in 1:N
    deg_a[a] = length(neighbors[a])
end

ϕ_links = zeros(Float64, (M, K+2))
ϕ_nonlinks=zeros(Float64, (N, K))
for m in 1:M
    ϕ_links[m, K+1] = links[m,1]
    ϕ_links[m, K+2] = links[m,2]
    for k in 1:K
        ϕ_links[m,k] = 1.0/K
    end
end
##watch out, data frame is more convenient as links are not Float64
# ϕ_nonlinks = update_ϕ_nonlinks(ϕ_links)
# for a in 1:N
#     for k in 1:K
#         ϕ_nonlinks[a,k]=(sum(ϕ_links[ϕ_links[:,K+1].== a,k]) + sum(ϕ_links[ϕ_links[:,K+2].== a,k]))/sum(adj_matrix[a,:])
#     end
# end
for k in 1:K, a in 1:N
    ϕ_nonlinks[a,k]=(sum(view(ϕ_links, view(ϕ_links, :,(K+1)) .== a, k)) + sum(view(ϕ_links, view(ϕ_links, :,(K+2)) .== a, k))) /deg_a[a]
end

# for k in 1:K, a in 1:N
#     ϕ_nonlinks[ϕ_links[a, K+1], k] += ϕ_links[a, k]
#     ϕ_nonlinks[ϕ_links[a, K+2], k] += ϕ_links[a, k]
# end
srand(1234)
γ = zeros(Float64, (N, K))
for k in 1:K, a in 1:N
    dunif = Uniform(0,1)
    γ[a,k] = rand(dunif)
end

total_pairs = convert(Int64, N*(N-1)/2)
ones_prob = M/total_pairs
τ0 = zeros(Float64, K)
τ1 = zeros(Float64, K)
for k in 1:K
    τ0[k] = η0
    τ1[k] = η1
end



MAX_ITER=1000
min_threshold = 1e-8


Elog_Θ = Elog_Dirichlet(γ)
Elog_β = Elog_Beta(τ0, τ1)


log_ϵ = log(ϵ)
log_1_minus_ϵ = log1p(-ϵ)
log_phi = zeros(Float64, (M,K))
CONVERGED= false
iteration=1
ELBO=Float64[]
"""
compute_elbo(log_ϵ, log_1_minus_ϵ, ϕ_links, ϕ_nonlinks, Elog_Θ, Elog_β, non_links, η0, η1, τ0, τ1, γ, α)
computes the variational lower bound(ELBO) for the E-step
# Arguments
* `log_ϵ`: is log of ϵ value
* `log_1_minus_ϵ`: is log of (1-ϵ) value
* `ϕ_links`: is a matrix of N×K variational parameters
* `ϕ_nonlinks`: is a matrix of N×K variational parameters
* `Elog_Θ`: is a matrix of N×K
* `Elog_β`: is a matrix of K×2
* `non_links`: the list of neighbors of each node
* `η0`: model scale hyperparameter of the Beta dsitributed parameter
* `η1`: model shape hyperparameter of the Beta dsitributed parameter
* `τ0`: scale hyperparameter of the Beta dsitributed variational parameter
* `τ1`: shape hyperparameter of the Beta dsitributed variational parameter
* `γ` : variational parameter for Dirichlet distributd Θ, a matrix of N×K
* `α` : model Dirichlet parameter for community strength
"""
function compute_elbo(log_ϵ, log_1_minus_ϵ, ϕ_links, ϕ_nonlinks, Elog_Θ, Elog_β, non_links, η0, η1, τ0, τ1, γ, α)
    num_nonlinks = size(non_links)[1]
    s = 0.0
    ##inefficient(define the view(A, ind) for dataframes, or use matrix form of links)
    for k in 1:K, m in 1:M
        s = s + ϕ_links[m,k] * (Elog_Θ[links[m,1],k]+Elog_Θ[links[m,2],k]+Elog_β[k,1]-log(ϕ_links[m,k]-log_ϵ))+log_ϵ
    end
    for k in 1:K, mn in 1:num_nonlinks
        x1 = ϕ_nonlinks[non_links[mn,1],k]
        x2 = ϕ_nonlinks[non_links[mn,2],k]
        s = s + x1*x2*(Elog_β[k,2]-log_1_minus_ϵ)+x1*(Elog_Θ[non_links[mn, 1],k]-log(x1))+x2*(Elog_Θ[non_links[mn, 2],k]-log(x2))+log_1_minus_ϵ
    end
    for a in 1:N
        s = s - lgamma(sum(γ[a,:]))
    end

    for k in 1:K
        s = s+(η0-τ0[k])*Elog_β[k,1]+(η1-τ1[k])*Elog_β[k,1]-lgamma(τ0[k]+τ1[k])+lgamma(τ0[k])+lgamma(τ1[k])
    end
    for k in 1:K, a in 1:N
         s = s + (α[k]-γ[a,k])*Elog_Θ[a,k]+lgamma(γ[a,k])
    end
    s
end


while !CONVERGED || (iteration > MAX_ITER)
    for m in 1:M
        sum_log_philink = 0
        for k in 1:K
            log_phi[m,k] = Elog_Θ[links[m,1],k]+Elog_Θ[links[m,2],k]+Elog_β[k,1]
            if k > 1
                sum_log_philink = log_sum(sum_log_philink, log_phi[m,k])
            else
                sum_log_philink = log_phi[m,k]
            end
        end
        for k in 1:K
            ϕ_links[m,k] = exp(log_phi[m,k] - sum_log_philink)
        end
    end
    #test
    for k in 1:K, a in 1:N
        ϕ_nonlinks[a,k]=(sum(ϕ_links[ϕ_links[:,K+1].== a,k]) + sum(ϕ_links[ϕ_links[:,K+2].== a,k]))/sum(adj_matrix[a,:])
    end
    # ϕ_nonlinks = update_ϕ_nonlinks(ϕ_links)
    #test
    for k in 1:K, a in 1:N
        γ[a,k] = α[k]+deg_a[a]*ϕ_nonlinks[a,k]+(N-1-deg_a[a])*ϕ_nonlinks[a,k]
    end
    for k in 1:K
        τ0[k] = η0 + sum(ϕ_links[:,k])
    end
    s1 = s2 = zeros(Float64, K)
    for k in 1:K
        s1[k] = sum(ϕ_nonlinks[:,k])
        for  a in 1:N
            s2[k] = s2[k] + ϕ_nonlinks[a,k]*ϕ_nonlinks[a,k]
        end
    end
    s3 = zeros(Float64, K)
    for m in 1:M
        for k in 1:K
            s3[k] = s3[k] + ϕ_nonlinks[links[m,1],k]+ ϕ_nonlinks[links[m,2],k]
        end
    end
    for k in 1:K
        τ1[k] = η1 + (s1[k]*s1[k] - s2[k])/2 - s3[k]
    end
    Elog_Θ = Elog_Dirichlet(γ)
    Elog_β = Elog_Beta(τ0, τ1)
    if rem(iteration, 1) == 0
        x = compute_elbo(log_ϵ, log_1_minus_ϵ, ϕ_links, ϕ_nonlinks, Elog_Θ,Elog_β, non_links, η0, η1, τ0, τ1,γ, α)
        print(iteration)
        print(": ")
        print(x)
        print("\n")
        push!(ELBO, x)
        if length(ELBO) > 1
            if abs(ELBO[length(ELBO)] - ELBO[length(ELBO)-1]) < min_threshold
                CONVERGED=true
                print("CONVERGED\n")
                break
            end
        end
    end
    iteration = iteration+1
    if iteration > MAX_ITER
        print("reached over MAX_ITER\n")
        break
    end
end


elbo_df = DataFrame(ITERATIONS=1:length(ELBO), ELBO_VALUE=ELBO[1:length(ELBO)])

ELBO
writetable("elbo_julia.csv", elbo_df)
elbo_df
Gadfly.plot(elbo_df, x="ITERATIONS", y="ELBO_VALUE",Geom.point,Geom.line)

end
