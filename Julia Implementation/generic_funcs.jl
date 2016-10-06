__precompile__
module GenericFuncs
export log_sum, get_static_neighbors, get_static_nonneighbors, get_links,get_nonlinks
export get_neighbors, get_nonneighbors, sort_by_argmax,Elog_Dirichlet,Elog_Beta
using DataFrames
function sort_by_argmax(mat)
    N=size(mat)[1]
    K=size(mat)[2]
    amax=zeros(Int64, N)
    for a in 1:N
        amax[a] = indmax(mat[a,:])
    end

    mat_tmp = similar(mat)
    count = 1
    for j in 1:maximum(amax)
        for i in 1:N
            if amax[i]==j
                mat_tmp[count,:]=mat[i,:]
                count =count+1
            end
        end
    end
    mat_tmp
end

function get_neighbors(mat, a)
    N = size(mat)[1]
    neigh_a = Int64[]
    for b in 1:N
        if mat[a,b] == 1
            push!(neigh_a,b)
        end
    end
    neigh_a
end

function get_nonneighbors(mat, a)
    N = size(mat)[1]
    neigh_a = get_neighbors(mat, a)
    nonneigh_a = Int64[]
    for i in 1:N
        count = 0
        for j in neigh_a
            if i == j
                count = count + 1
            end
        end
        if count == 0
            push!(nonneigh_a,i)
        end
    end
    nonneigh_a
end
function get_links(mat)
    N = size(mat)[1]
    num_links = convert(Int64, sum(mat)/2)
    df = DataFrame(X1=1:num_links, X2=1:num_links)
    l_count=1
    for a in 1:N
        for b in a:N
            if b == a
                continue
            end
            if mat[a,b] == 1
                df[l_count,1] = a
                df[l_count,2] = b
                l_count = l_count + 1
            end
        end
    end
    df
end

function get_nonlinks(mat)
    N=size(mat)[1]
    num_links = convert(Int64, sum(mat)/2)
    num_nonlinks = convert(Int64, N*(N-1)/2 - num_links)
    df = DataFrame(X1=1:num_nonlinks, X2=1:num_nonlinks)
    nl_count=1
    for a in 1:N
        for b in a:N
            if b == a
                continue
            end
            if mat[a,b] == 0
                df[nl_count,1] = a
                df[nl_count,2] = b
                nl_count = nl_count + 1
            end
        end
    end
    df
end
function get_static_neighbors(mat)
    N = size(mat)[1]
    stat_neigh = Array(Vector{Int64}, N)
    for a in 1:N
        stat_neigh[a] = get_neighbors(mat, a)
    end
    stat_neigh
end


function get_static_nonneighbors(mat)
    N= size(mat)[1]
    stat_nonneigh = Array(Vector{Int64}, N)
    for a in 1:N
        stat_nonneigh[a] = get_nonneighbors(mat, a)
    end
    stat_nonneigh
end

function log_sum(log_a, log_b)
    v = 0.0
    if log_a < log_b
        if exp(log_a - log_b) < 1e-10
            v = log_b + log1p(exp(log_a-log_b))
        else
            v = log_b + log(1 + exp(log_a-log_b))
        end
    else
        if exp(log_b - log_a) < 1e-10
            v = log_a + log1p(exp(log_b - log_a))
        else
            v = log_a + log(1 + exp(log_b - log_a))
        end
    end
    v
end
function Elog_Dirichlet(mat)
    N = size(mat)[1]
    K = size(mat)[2]
    temp = zeros(Float64, (N,K))
    for a in 1:N
        for k in 1:K
            temp[a,k] = digamma(mat[a,k]) - digamma(sum(mat[a,:]))
        end
    end
    temp
end
function Elog_Beta(l0, l1)
    K = length(l1)
    temp = zeros(Float64, (K,2))
    for k in 1:K
        temp[k,1] = digamma(l0[k]) - digamma(l0[k]+l1[k])
        temp[k,2] = digamma(l1[k]) - digamma(l0[k]+l1[k])
    end
    temp
end
end
