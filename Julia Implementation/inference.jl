__precompile__
module Inference
export ELBO
# include("generic_funcs.jl")
# include("genr_network.jl")
# include("work__space.jl")
# include("local_update.jl")
using DataFrames
using Distributions
using Gadfly
using GenericFuncs:get_links, get_nonlinks, get_static_neighbors, get_static_nonneighbors,Elog_Dirichlet,Elog_Beta, log_sum
using GenrNetwork
using LocalUpdate:update_ϕ_nonlinks


M = convert(Int64, sum(adj_matrix/2))
links=get_links(adj_matrix)
non_links = get_nonlinks(adj_matrix)
neighbors = get_static_neighbors(adj_matrix)
non_neighbors = get_static_nonneighbors(adj_matrix)

ϕ_links = zeros(Float64, (M, K+2))
for m in 1:M
    ϕ_links[m, K+1] = links[m,1]
    ϕ_links[m, K+2] = links[m,2]
    for k in 1:K
        ϕ_links[m,k] = 1.0/K
    end
end
##watch out, data frame is more convenient as links are not Float64
ϕ_nonlinks = update_ϕ_nonlinks(ϕ_links)

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

deg_a = zeros(Int64, N)
for a in 1:N
    deg_a[a] = length(neighbors[a])
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

function compute_elbo(log_ϵ, log_1_minus_ϵ, ϕ_links, ϕ_nonlinks, Elog_Θ, Elog_β, non_links, η0, η1, τ0, τ1, γ, α)
    num_nonlinks = size(non_links)[1]
    s = 0.0
    for m in 1:M
        for k in 1:K
            s = s + ϕ_links[m,k] * (Elog_Θ[links[m,1],k]+Elog_Θ[links[m,2],k]+Elog_β[k,1]-log(ϕ_links[m,k]-log_ϵ))+log_ϵ
        end
    end
    for mn in 1:num_nonlinks
        for k in 1:K
            x1 = ϕ_nonlinks[non_links[mn,1],k]
            x2 = ϕ_nonlinks[non_links[mn,2],k]
            s = s + x1*x2*(Elog_β[k,2]-log_1_minus_ϵ)+x1*(Elog_Θ[non_links[mn, 1],k]-log(x1))+x2*(Elog_Θ[non_links[mn, 2],k]-log(x2))+log_1_minus_ϵ
        end
    end
    for a in 1:N
        s = s - lgamma(sum(γ[a,:]))
    end

    for k in 1:K
            s = s+(η0-τ0[k])*Elog_β[k,1]+(η1-τ1[k])*Elog_β[k,1]-lgamma(τ0[k]+τ1[k])+lgamma(τ0[k])+lgamma(τ1[k])
    end
    for a in 1:N
        for k in 1:K
             s = s + (α[k]-γ[a,k])*Elog_Θ[a,k]+lgamma(γ[a,k])
        end
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
    ϕ_nonlinks = update_ϕ_nonlinks(ϕ_links)
    #test
    for a in 1:N
        for k in 1:K
            γ[a,k] = α[k]+deg_a[a]*ϕ_nonlinks[a,k]+(N-1-deg_a[a])*ϕ_nonlinks[a,k]
        end
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
elbo_df
p= Gadfly.plot(elbo_df, x="ITERATIONS", y="ELBO_VALUE",Geom.point,Geom.line)
end
