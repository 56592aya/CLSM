using DataFrames
using Distributions
using LightGraphs
#using Plots
#using StatPlots
#using Vega
#using Gadfly
#using PyPlot
srand(1234)
N=100;K=4;ϵ = 1e-30;η0 = 10;η1 = 1;
##Creating β################################
β = zeros(Float64, (K,K))
dbeta=Beta(η0, η1)
β_diag = rand(dbeta, K)
for k in 1:K
    β[k,k] = β_diag[k]
end
β=β + ϵ*(ones(Int64, (K,K))-eye(Int64, K))

###########################################
#Creating α################################
α=Array(Float64, K)
for k in 1:K
    α[k] = 0.01
end
ddir = Dirichlet(α)
Θ = zeros(Float64, (N,K))
for a in 1:N
    Θ[a,:] = rand(ddir)
end

function sort_by_argmax!(mat::Array{Float64, 2})
    n_row=size(mat)[1]
    ind_max=zeros(Int64, n_row)
    for a in 1:n_row
        ind_max[a] = indmax(mat[a,:])
    end
    mat_tmp = similar(mat)
    count = 1
    for j in 1:maximum(ind_max)
        for i in 1:n_row
            if ind_max[i] == j
                mat_tmp[count,:] = mat[i,:]
                count = count + 1
            end
        end
    end
    mat[:]=mat_tmp[:] ###This is important in arrays
    mat
end

sort_by_argmax!(Θ)
#################################################
###Creating the adjacency matrix#################
adj_matrix = zeros(Int64, (N,N))

for b in 1:N
    for a in 1:N
        if b != a
            dmult_a = Multinomial(1,vec(Θ[a,:]))
            dmult_b = Multinomial(1,vec(Θ[b,:]))
            z_send_ab = rand(dmult_a,1)
            z_recv_ab = rand(dmult_b,1)
            z_send_ab_idx = indmax(z_send_ab)
            z_recv_ab_idx = indmax(z_recv_ab)
            dbinom = Binomial(1,β[z_send_ab_idx, z_recv_ab_idx])
            adj_matrix[a,b] = rand(dbinom, 1)[1]
        end
    end
end
##Creating the DiGraph
network  = DiGraph(N)
for b in 1:N
    for a in 1:N
        if adj_matrix[a,b] == 1
            add_edge!(network, a, b)
        end
    end
end

links = edges(network)
M_total = network.ne
sinks = fadj(network)
sources = badj(network)
###gotta remove self index

non_sinks=[setdiff(deleteat!(Vector(1:N), index), value) for (index,value) in enumerate(sinks) ]
non_sources=[setdiff(deleteat!(Vector(1:N), index), value) for (index,value) in enumerate(sources) ]
out_degrees = map(x->outdegree(network,x), vertices(network))
in_degrees = map(x->indegree(network,x), vertices(network))


##################################################
####Initializing variational params###############
#ϕ_links_send[a,b] means link is from a to b
#ϕ_links_recv[a,b] means link is from a to b
ϕ_links_send=DataFrame()
for i in 1:K
    ϕ_links_send[Symbol("K$i")]=repeat([1.0/K], outer = [M_total])
end
ϕ_links_send[:e]=[v for (i,v) in enumerate(edges(network))]
ϕ_links_send[:sourceAt]=[v.first for (i,v) in enumerate(edges(network))]
ϕ_links_send[:sinkAt]=[v.second for (i,v) in enumerate(edges(network))]
#####

ϕ_links_recv=DataFrame()
for i in 1:K
    ϕ_links_recv[Symbol("K$i")]=repeat([1.0/K], outer = [M_total])
end
ϕ_links_recv[:e]=[v for (i,v) in enumerate(edges(network))]
ϕ_links_recv[:sourceAt]=[v.first for (i,v) in enumerate(edges(network))]
ϕ_links_recv[:sinkAt]=[v.second for (i,v) in enumerate(edges(network))]
###
ϕ_nonlinks_send=DataFrame()
for i in 1:K
    ϕ_nonlinks_send[Symbol("K$i")]=repeat([1.0/K], outer = [N])
end
ϕ_nonlinks_recv=DataFrame()
for i in 1:K
    ϕ_nonlinks_recv[Symbol("K$i")]=repeat([1.0/K], outer = [N])
end
#################
γ = zeros(Float64, (N,K))
for k in 1:K
    for a in 1:N
        #γ[a,k] = rand()
        γ[a,k] = 2.0*N/K
    end
end

#############################
τ=zeros(2,K)
τ[1,:] = repeat([η0], outer=[K])
τ[2,:] = repeat([η1], outer=[K])
#############################
function log_sum(log_a::Float64, log_b::Float64)
    out = 0.0
    if log_a < log_b
        if exp(log_a - log_b) < 1e-10
            out = log_b + log1p(exp(log_a-log_b))
        else
            out = log_b + log(1 + exp(log_a-log_b))
        end
    else
        if exp(log_b - log_a) < 1e-10
            out = log_a + log1p(exp(log_b - log_a))
        else
            out = log_a + log(1 + exp(log_b - log_a))
        end
    end
    out
end
########################
function Elog_Dirichlet(mat::Array{Float64, 2})
    N = size(mat)[1]
    K = size(mat)[2]
    out = similar(mat)
    for a in 1:N
        for k in 1:K
            out[a,k] = digamma(mat[a,k]) - digamma(sum(mat[a,:]))
        end
    end
    out
end
#########################
function Elog_Beta(mat::Array{Float64, 2})
    K = size(mat)[2]
    out = similar(mat)
    for k in 1:K
        out[1,k] = digamma(mat[1,k]) - digamma(sum(mat[:,k]))
        out[2,k] = digamma(mat[2,k]) - digamma(sum(mat[:,k]))
    end
    out
end

###################
Elog_Θ = zeros(Float64, (N,K))
Elog_Θ[:] = Elog_Dirichlet(γ)
Elog_β=zeros(Float64, (2,K))
Elog_β[:] = Elog_Beta(τ)

##################################

all_pairs = Array{Pair{Int64, Int64},1}(N*N-N)
count = 1
for a in 1:N
    for b in 1:N
        if b != a
            all_pairs[count] = a=>b
            count +=1
        end
    end
end
nonlinks=setdiff(all_pairs, ϕ_links_send[:e])
function compute_elbo(links, ϕ_links_send,ϕ_links_recv,Elog_β,log_ϵ,Elog_Θ,η0,η1,τ,α,γ,nonlinks, log_1_minus_ϵ,ϕ_nonlinks_send, ϕ_nonlinks_recv)
    s = 0
    for (index,value) in enumerate(links)
        for k in 1:K
            s = s+(ϕ_links_send[index, Symbol("K$k")]*ϕ_links_recv[index, Symbol("K$k")])*(Elog_β[1,k]-log_ϵ)+log_ϵ+
            ϕ_links_send[index, Symbol("K$k")]*(Elog_Θ[value.first,k]-log(ϕ_links_send[index, Symbol("K$k")]))+
            ϕ_links_recv[index, Symbol("K$k")]*(Elog_Θ[value.second,k]-log(ϕ_links_recv[index, Symbol("K$k")]))
        end
    end
    for k in 1:K
        s+=(η0-τ[1,k])*Elog_β[1,k] + (η1-τ[2,k])*Elog_β[2,k]-lgamma(sum(τ[:,k]))+lgamma(τ[1,k])+lgamma(τ[2,k])
    end
    for k in 1:K
        for a in 1:N
            s += (α[k]-γ[a,k])*Elog_Θ[a,k]+lgamma(γ[a,k])
        end
    end
    for a in 1:N
        s -= lgamma(sum(γ[a,:]))
    end

    for m in 1:length(nonlinks)
        for k in 1:K
            s+=(ϕ_nonlinks_send[nonlinks[m].first,Symbol("K$k")]*ϕ_nonlinks_recv[nonlinks[m].second,Symbol("K$k")])*(Elog_β[2,k]-log_1_minus_ϵ)+log_1_minus_ϵ+
            ϕ_nonlinks_send[nonlinks[m].first,Symbol("K$k")]*(Elog_Θ[nonlinks[m].first, k]-ϕ_nonlinks_send[nonlinks[m].first,Symbol("K$k")]) +
            ϕ_nonlinks_recv[nonlinks[m].second,Symbol("K$k")]*(Elog_Θ[nonlinks[m].second, k]-ϕ_nonlinks_recv[nonlinks[m].second,Symbol("K$k")])
        end
    end
    s
end
#compute_elbo(links, ϕ_links_send,ϕ_links_recv,Elog_β,log_ϵ,Elog_Θ,η0,η1,τ,α,γ,nonlinks, log_1_minus_ϵ,ϕ_nonlinks_send, ϕ_nonlinks_recv)
MAX_ITER=1000###FOR NOW
min_threshold = 1e-8
log_ϵ = log(ϵ)
log_1_minus_ϵ = log1p(-ϵ)
CONVERGED= false
iteration=1
ELBO=Float64[]
log_phi_send = zeros(Float64, (M_total,K))
log_phi_recv = zeros(Float64, (M_total,K))
s0 = zeros(Float64, K)
s1 = zeros(Float64, K)
s2 = zeros(Float64, K)
s3 = zeros(Float64, (N,K))

while !CONVERGED || (iteration > MAX_ITER)
    for (index,value) in enumerate(links)
        sum_log_philink_send = 0
        sum_log_philink_recv = 0
        for k in 1:K
            log_phi_send[index,k] = ϕ_links_recv[index,k]*(Elog_β[1,k]-log_ϵ)+Elog_Θ[value.first,k]
            log_phi_recv[index,k] = ϕ_links_send[index,k]*(Elog_β[1,k]-log_ϵ)+Elog_Θ[value.second,k]
            if k > 1
                sum_log_philink_send = log_sum(sum_log_philink_send, log_phi_send[index,k])
                sum_log_philink_recv = log_sum(sum_log_philink_recv, log_phi_recv[index,k])
            else
                sum_log_philink_send = log_phi_send[index,k]
                sum_log_philink_recv = log_phi_recv[index,k]
            end
        end
        for k in 1:K
            ϕ_links_send[index,k] = exp(log_phi_send[index,k] - sum_log_philink_send)
            ϕ_links_recv[index,k] = exp(log_phi_recv[index,k] - sum_log_philink_recv)
        end
    end
    writetable("ϕ_links_send_$iteration.csv",ϕ_links_send)
    writetable("ϕ_links_recv_$iteration.csv",ϕ_links_recv)
    for k in 1:K
        for a in 1:N
            #a as non-source
            ϕ_nonlinks_send[a,Symbol("K$k")] = (sum(ϕ_links_send[ϕ_links_send[:sourceAt] .== a,Symbol("K$k")]) + sum(ϕ_links_recv[ϕ_links_recv[:sinkAt] .== a,Symbol("K$k")]))/(in_degrees[a]+out_degrees[a])
            #a as non-sink
            ϕ_nonlinks_recv[a,Symbol("K$k")] = ϕ_nonlinks_send[a,Symbol("K$k")]
        end
    end
    writetable("ϕ_nonlinks_send_$iteration.csv",ϕ_nonlinks_send)
    writetable("ϕ_nonlinks_recv_$iteration.csv",ϕ_nonlinks_recv)
    for k in 1:K
        for a in 1:N
            γ[a,k] = α[k]+(sum(ϕ_links_send[ϕ_links_send[:sourceAt] .== a,Symbol("K$k")])+sum(ϕ_links_recv[ϕ_links_recv[:sinkAt] .== a,Symbol("K$k")]))*(2*(N-1))/(in_degrees[a]+out_degrees[a])
        end
    end
    writedlm("gamma_$iteration.csv", γ)
    s0 = zeros(Float64, K)
    for k in 1:K
        for m in 1:M_total
            s0[k] += ϕ_links_send[m,Symbol("K$k")]*ϕ_links_recv[m,Symbol("K$k")]
        end
    end

    s1 = zeros(Float64, K)
    for k in 1:K
         for a in 1:N
             s1[k]+=ϕ_nonlinks_send[a,Symbol("K$k")]
         end
    end
    s3 = zeros(Float64, (N,K))
    for a in 1:N
        for k in 1:K
            for b in sinks[a]
                s3[a,k]+=ϕ_nonlinks_recv[b,Symbol("K$k")]
            end
        end
    end

    s2 = zeros(Float64, K)
    for a in 1:N
        for k in 1:K
            s2[k] += ϕ_nonlinks_send[a,Symbol("K$k")]*(ϕ_nonlinks_send[a,Symbol("K$k")]+s3[a,k])
        end
    end

    for k in 1:K
        τ[1,k] = η0 + s0[k]
        τ[2,k] = η1 + s1[k]*s1[k] -s2[k]
    end
    writedlm("tau_$iteration.csv", τ)
    Elog_Θ[:] = Elog_Dirichlet(γ)
    Elog_β[:] = Elog_Beta(τ)
    writedlm("ELOG_THETA_$iteration.csv",Elog_Θ)
    writedlm("ELOG_BETA_$iteration.csv",Elog_β)
    if rem(iteration, 10) == 0
        x = compute_elbo(links, ϕ_links_send,ϕ_links_recv,Elog_β,log_ϵ,Elog_Θ,η0,η1,τ,α,γ,nonlinks, log_1_minus_ϵ,ϕ_nonlinks_send, ϕ_nonlinks_recv)
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
