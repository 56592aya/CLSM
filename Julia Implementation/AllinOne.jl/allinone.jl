module AllInOne
using DataFrames
using Distributions
using Gadfly
using LightGraphs
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
###Network
N=200
K=8
ϵ = 1e-30
η0 = 10
η1 = 1
srand(1234)
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


network  = Graph(N)
for a in 1:N
    for b in a:N
        if adj_matrix[a,b] == 1
            add_edge!(network, a, b)
        end
    end
end

links = edges(network)
M = network.ne

sinks = fadj(network)
sources = badj(network)

non_links = get_nonlinks(adj_matrix)
non_neighbors = get_static_nonneighbors(adj_matrix)
deg_a = map(length, neighbors)

ϕ_links = convert(DataFrame, zeros(Float64, (M, K+2)))
ϕ_links[:,K+1] = convert(Array{Int64,1}, ϕ_links[:,K+1])
ϕ_links[:,K+2] = convert(Array{Int64,1}, ϕ_links[:,K+2])

ϕ_nonlinks=zeros(Float64, (N, K))
#links.m

for (index, value) in enumerate(links)
    ϕ_links[index, K+1] = value.first
    ϕ_links[index, K+2] = value.second
    for k in 1:K
        ϕ_links[index,k] = 1.0/K
    end
end


for k in 1:K, m in 1:M
    ϕ_nonlinks[ϕ_links[m,K+1],k] += ϕ_links[m,k]
    ϕ_nonlinks[ϕ_links[m,K+2],k] += ϕ_links[m,k]
end
for k in 1:K, a in 1:N
    ϕ_nonlinks[a,k] = ϕ_nonlinks[a,k]/deg_a[a]

srand(1234)
γ = zeros(Float64, (N, K))
for k in 1:K, a in 1:N
    dunif = Uniform(0,1)
    γ[a,k] = rand(dunif)
end

total_pairs = binomial(N,2)

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
    for k in 1:K, (index, value) in enumerate(links)
        s = s + ϕ_links[index,k] * (Elog_Θ[value.first,k]+Elog_Θ[value.second,k]+Elog_β[k,1]-log(ϕ_links[index,k]-log_ϵ))+log_ϵ
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
    for (index, value) in enumerate(links)
        sum_log_philink = 0
        for k in 1:K
            log_phi[index,k] = Elog_Θ[value.first,k]+Elog_Θ[value.second,k]+Elog_β[k,1]
            if k > 1
                sum_log_philink = log_sum(sum_log_philink, log_phi[index,k])
            else
                sum_log_philink = log_phi[index,k]
            end
        end
        for k in 1:K
            ϕ_links[index,k] = exp(log_phi[index,k] - sum_log_philink)
        end
    end
    #test
    #for k in 1:K, a in 1:N
    #    ϕ_nonlinks[a,k]=(sum(ϕ_links[ϕ_links[:,K+1].== a,k]) + sum(ϕ_links[ϕ_links[:,K+2].== a,k]))/sum(adj_matrix[a,:])
    #end
    for k in 1:K, m in 1:M
        ϕ_nonlinks[ϕ_links[m,K+1],k] += ϕ_links[m,k]
        ϕ_nonlinks[ϕ_links[m,K+2],k] += ϕ_links[m,k]
    end
    for k in 1:K, a in 1:N
        ϕ_nonlinks[a,k] = ϕ_nonlinks[a,k]/deg_a[a]
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
    for (index, value) in enumerate(links)
        for k in 1:K
            s3[k] = s3[k] + ϕ_nonlinks[value.first,k]+ ϕ_nonlinks[value.second,k]
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

estimate_beta(τ0, τ1)
diag(β)





end
