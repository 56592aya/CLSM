__precompile__
module GenericFuncs

export log_sum, get_static_neighbors, get_static_nonneighbors, get_links,
get_nonlinks, get_neighbors, get_nonneighbors, sort_by_argmax,
Elog_Dirichlet,Elog_Beta

using DataFrames

"""
sort_by_argmax(mat)
sorts a matrix `mat` by their argmax positions along the rows.
# Arguments
* `mat`: a two dimensional array.
# Examples

```jldoctest
julia> a = [1 2;4 0;2 3]
julia> sort_by_argmax(a)
```
"""
function sort_by_argmax(mat::Array{Float64, 2})
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

"""
get_neighbors(mat,a)
finds the neighbors(1 positions) of `a` in `mat`
# Arguments
* `mat`: a square matrix.
* `a` : a single node
# Examples

```jldoctest
julia> a = 2
julia> mat = [0 1 1 0;1 0 0 0; 1 0 0 1;0 0 1 0]
julia> get_neighbors(mat, a)
```
"""
function get_neighbors(mat::Array{Int64, 2}, a)
    N = size(mat)[1]
    neigh_a = Int64[]
    for b in 1:N
        if mat[a,b] == 1
            push!(neigh_a,b)
        end
    end
    neigh_a
end

"""
get_nonneighbors(mat,a)
finds the nonneighbors(0 positions) of `a` in `mat`
# Arguments
* `mat`: a square matrix.
* `a` : a single node
# Examples

```jldoctest
julia> a = 2
julia> mat = [0 1 1 0;1 0 0 0; 1 0 0 1;0 0 1 0]
julia> get_nonneighbors(mat, a)
```
"""
function get_nonneighbors(mat::Array{Int64, 2}, a::Int64)
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

"""
get_links(mat)
gets all the unique links in a square matrix `mat`
# Arguments
* `mat`: a square matrix.
# Examples

```jldoctest
julia> mat = [0 1 1 0;1 0 0 0; 1 0 0 1;0 0 1 0]
julia> get_links(mat)
```
"""
function get_links(mat::Array{Int64, 2})
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

"""
get_nonlinks(mat)
gets all the unique nonlinks in a square matrix `mat`
# Arguments
* `mat`: a square matrix.
# Examples

```jldoctest
julia> mat = [0 1 1 0;1 0 0 0; 1 0 0 1;0 0 1 0]
julia> get_nonlinks(mat)
```
"""
function get_nonlinks(mat::Array{Int64, 2})
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

"""
get_static_neighbors(mat)
creates all the neighbors lists for all nodes.
# Arguments
* `mat`: a square matrix.
# Examples

```jldoctest
julia> mat = [0 1 1 0;1 0 0 0; 1 0 0 1;0 0 1 0]
julia> get_static_neighbors(mat)
```
"""
function get_static_neighbors(mat::Array{Int64, 2})
    N = size(mat)[1]
    stat_neigh = Array(Vector{Int64}, N)
    for a in 1:N
        stat_neigh[a] = get_neighbors(mat, a)
    end
    stat_neigh
end

"""
get_static_neighbors(mat)
creates all the nonneighbors lists for all nodes.
# Arguments
* `mat`: a square matrix.
# Examples

```jldoctest
julia> mat = [0 1 1 0;1 0 0 0; 1 0 0 1;0 0 1 0]
julia> get_static_nonneighbors(mat)
```
"""
function get_static_nonneighbors(mat::Array{Int64, 2})
    N= size(mat)[1]
    stat_nonneigh = Array(Vector{Int64}, N)
    for a in 1:N
        stat_nonneigh[a] = get_nonneighbors(mat, a)
    end
    stat_nonneigh
end

"""
log_sum(log_a, log_b)
computes the `log(a+b)` given `log(a)` and `log(b)`.
# Arguments
* `log_a`: a `Float64` number, `log(a)`
* `log_b`: a `Float64` number, `log(b)`
# Examples

```jldoctest
julia> a = 10;b=5;
julia> log_a=log(a); log_b= log(b)
julia> log_sum(log_a, log_b)
```
"""
# arr = zeros(Int64, (10,2))
# out = similar(arr)
# x   = similar(arr)
# arr[:,1] = [1 2 3 4 5 6 7 8 9 10]
# arr[:,2] = [2 3 2 3 4 4 5 6 7 5]
# for j in 1:2, i in 1:10
#     x[i,j]=sum(view(arr,view(arr,:,2) .== i, j))
# end
# x
# x2 = zeros(Int, 10, 2)
# for j in 1:size(arr, 2), i in 1:size(arr, 1)
#     x2[arr[i, 2], j] += arr[i, j]
# end
# x2

function log_sum(log_a::Float64, log_b::Float64)
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

"""
Elog_Dirichlet(mat)
computes the expectation of log of a Dirichlet distributed parameter.
this is for each element at [a,k]: a={1,...,N}, K={1,...,K}, Ψ(γₐₖ)-Ψ(∑ₖγₐ,)
# Arguments
* `mat` : a two dimensional array of N×K
# Examples

```jldoctest
julia> mat=[0.5 0.6; 0.2 0.8; 0.9 0.1]
julia> Elog_Dirichlet(mat)
```
"""
function Elog_Dirichlet(mat::Array{Float64, 2})
    N = size(mat)[1]
    K = size(mat)[2]
    temp = zeros(Float64, (N,K))
    for a in 1:N
        for k in 1:K
            # Ψ(γₐₖ)-Ψ(∑ₖγₐ,)
            temp[a,k] = digamma(mat[a,k]) - digamma(sum(mat[a,:]))
        end
    end
    temp
end

"""
Elog_Beta`(mat)
computes the expectation of log of a Beta distributed parameter.
this is for each element at [k,i]: k={1,...,K}, i={1,2}, Ψ(τₖᵢ)-Ψ(τₖ₁+τₖ₂)
# Arguments
* `mat` : a two dimensional array of K×2
# Examples

```jldoctest
julia> mat=[0.5 0.6; 0.2 0.8; 0.9 0.1]
julia> Elog_Beta(mat)
```
"""
function Elog_Beta(l0::Array{Float64, 1}, l1::Array{Float64, 1})
    K = length(l1)
    temp = zeros(Float64, (K,2))
    for k in 1:K
        # Ψ(τₖᵢ)-Ψ(τₖ₁+τₖ₂)
        temp[k,1] = digamma(l0[k]) - digamma(l0[k]+l1[k])
        temp[k,2] = digamma(l1[k]) - digamma(l0[k]+l1[k])
    end
    temp
end

function estimate_beta(tau0::Array{Float64, 1}, tau1::Array{Float64, 1})
    K=length(tau0)
    beta_rate= zeros(Float64, K)
    for k in 1:K
        sum=0
        sum=sum+tau0[k]+tau1[k]
        beta_rate[k]=tau0[k]/sum
    end
    beta_rate
end

end
