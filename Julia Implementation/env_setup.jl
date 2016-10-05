module EnvSetup

using GenericFuncs:get_static_neighbors, get_links, get_static_nonneighbors, get_nonlinks, Elog_Dirichlet, Elog_Beta
importall GenrNetwork
using GlobalInit:γ, τ0, τ1



links=get_links(adj_matrix)
non_links = get_nonlinks(adj_matrix)
neighbors = get_static_neighbors(adj_matrix)
non_neighbors = get_static_nonneighbors(adj_matrix)

M=size(links)[1]
deg_a = zeros(Int64, N)

for a in 1:N
    deg_a[a] = length(neighbors[a])
end

MAX_ITER=1000
min_threshold = 1e-8


Elog_Θ = Elog_Dirichlet(γ)
Elog_β = Elog_Beta(τ0, τ1)



export Elog_Θ, Elog_β
export deg_a,MAX_ITER, min_threshold, M
export links, non_links, neighbors, non_neighbors


end
