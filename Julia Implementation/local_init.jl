module LocalInit
using GenrNetwork:K, N, adj_matrix, η0, η1
using Distributions
using DataFrames
using GenericFuncs:get_links
import LocalUpdate:update_ϕ_nonlinks

M = convert(Int64, sum(adj_matrix/2))
links=get_links(adj_matrix)
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
export ϕ_links, ϕ_nonlinks
end
