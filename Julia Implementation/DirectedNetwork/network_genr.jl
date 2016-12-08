__precompile__()

module NetworkGenr
using DataFrames
using Distributions
using LightGraphs

export N, K, ϵ, η0, η1, β, α,Θ, adj_matrix, network, links, links_total, sinks, sources,non_sinks,non_sources,in_degrees, out_degrees,all_pairs, nonlinks
export sort_by_argmax!

#################################################
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
#################################################
srand(1234)
N=200;K=4;ϵ = 1e-30;η0 = 10;η1 = 1;
##Creating β ################################
β = zeros(Float64, (K,K))
dbeta=Beta(η0, η1)
β_diag = rand(dbeta, K)
for k in 1:K
    β[k,k] = β_diag[k]
end
β=β + ϵ*(ones(Int64, (K,K))-eye(Int64, K))

###########################################
#Creating α ################################
α=Array(Float64, K)
for k in 1:K
    α[k] = 0.005
end
ddir = Dirichlet(α)
Θ = zeros(Float64, (N,K))
for a in 1:N
    Θ[a,:] = rand(ddir)
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
#################################################
##Creating the DiGraph###########################
network  = DiGraph(N)
for b in 1:N
    for a in 1:N
        if adj_matrix[a,b] == 1
            add_edge!(network, a, b)
        end
    end
end
#################################################
links = edges(network)
links_total = network.ne
sinks = fadj(network)
sources = badj(network)
#################################################
non_sinks=[setdiff(deleteat!(Vector(1:N), index), value) for (index,value) in enumerate(sinks) ]
non_sources=[setdiff(deleteat!(Vector(1:N), index), value) for (index,value) in enumerate(sources) ]
out_degrees = map(x->outdegree(network,x), vertices(network))
in_degrees = map(x->indegree(network,x), vertices(network))
#################################################
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
nonlinks=setdiff(all_pairs, links)
#################################################
#################################################

end
