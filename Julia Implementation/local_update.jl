module LocalUpdate
using GenrNetwork:K, N, adj_matrix

function update_ϕ_nonlinks(mat)
    temp = zeros(Float64, (N,K))
    for a in 1:N
        for k in 1:K
            temp[a,k]=(sum(mat[mat[:,K+1].== a,k]) + sum(mat[mat[:,K+2].== a,k]))/sum(adj_matrix[a,:])
        end
    end
    temp
end
export update_ϕ_nonlinks
end
