__precompile__()
module Youtilities
using DataFrames
using NetworkGenr

export Elog_Beta, Elog_Dirichlet, log_sum,compute_elbo

#################################################
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
#################################################
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
#################################################
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
#################################################
function Elog_Beta(mat::Array{Float64, 2})
    K = size(mat)[2]
    out = similar(mat)
    for k in 1:K
        out[1,k] = digamma(mat[1,k]) - digamma(sum(mat[:,k]))
        out[2,k] = digamma(mat[2,k]) - digamma(sum(mat[:,k]))
    end
    out
end
#################################################
#################################################

end
