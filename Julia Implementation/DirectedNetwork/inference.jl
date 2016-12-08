__precompile__()
module Inference
using DataFrames
using LightGraphs
using NetworkGenr
using Youtilities
export ϕ_links_send,ϕ_links_recv,ϕ_nonlinks_send,ϕ_nonlinks_recv,γ,τ,ELBO
######Initializations
#################################################
#################################################
ϕ_links_send=DataFrame()
for i in 1:K
    ϕ_links_send[Symbol("K$i")]=repeat([1.0/K], outer = [links_total])
end
ϕ_links_send[:e]=[v for (i,v) in enumerate(edges(network))]
ϕ_links_send[:sourceAt]=[v.first for (i,v) in enumerate(edges(network))]
ϕ_links_send[:sinkAt]=[v.second for (i,v) in enumerate(edges(network))]
#################################################
#################################################
ϕ_links_recv=DataFrame()
for i in 1:K
    ϕ_links_recv[Symbol("K$i")]=repeat([1.0/K], outer = [links_total])
end
ϕ_links_recv[:e]=[v for (i,v) in enumerate(edges(network))]
ϕ_links_recv[:sourceAt]=[v.first for (i,v) in enumerate(edges(network))]
ϕ_links_recv[:sinkAt]=[v.second for (i,v) in enumerate(edges(network))]
#################################################
ϕ_nonlinks_send=DataFrame()
for i in 1:K
    ϕ_nonlinks_send[Symbol("K$i")]=repeat([1.0/K], outer = [N])
end
#################################################
#################################################
ϕ_nonlinks_recv=DataFrame()
for i in 1:K
    ϕ_nonlinks_recv[Symbol("K$i")]=repeat([1.0/K], outer = [N])
end
#################################################
#################################################
γ = zeros(Float64, (N,K))
for k in 1:K
    for a in 1:N
        γ[a,k] = rand()
        #γ[a,k] = 2.0*N/K
    end
end
#################################################
#################################################
τ=zeros(2,K)
τ[1,:] = repeat([η0], outer=[K])
τ[2,:] = repeat([η1], outer=[K])
#################################################
#################################################
#################################################

Elog_Θ = zeros(Float64, (N,K))
Elog_Θ[:] = Elog_Dirichlet(γ)
Elog_β=zeros(Float64, (2,K))
Elog_β[:] = Elog_Beta(τ)
#################################################
#################################################

MAX_ITER=1000###FOR NOW
min_threshold = 1e-8
log_ϵ = log(ϵ)
log_1_minus_ϵ = log1p(-ϵ)
CONVERGED= false
iteration=1
ELBO=Float64[]
log_phi_send = zeros(Float64, (links_total,K))
log_phi_recv = zeros(Float64, (links_total,K))
s0 = zeros(Float64, K)
s1 = zeros(Float64, K)
s2 = zeros(Float64, K)
s3 = zeros(Float64, (N,K))
#################################################
#################################################

#################################################

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

    # writetable("TestPrints/ϕ_links_send_$iteration.csv",ϕ_links_send)
    # writetable("TestPrints/ϕ_links_recv_$iteration.csv",ϕ_links_recv)
    for k in 1:K
        for a in 1:N
            #a as non-source
            ϕ_nonlinks_send[a,Symbol("K$k")] = (sum(ϕ_links_send[ϕ_links_send[:sourceAt] .== a,Symbol("K$k")]) + sum(ϕ_links_recv[ϕ_links_recv[:sinkAt] .== a,Symbol("K$k")]))/(in_degrees[a]+out_degrees[a])
            #a as non-sink
            ϕ_nonlinks_recv[a,Symbol("K$k")] = ϕ_nonlinks_send[a,Symbol("K$k")]
        end
    end
    # writetable("TestPrints/ϕ_nonlinks_send_$iteration.csv",ϕ_nonlinks_send)
    # writetable("TestPrints/ϕ_nonlinks_recv_$iteration.csv",ϕ_nonlinks_recv)
    for k in 1:K
        for a in 1:N
            γ[a,k] = α[k]+(sum(ϕ_links_send[ϕ_links_send[:sourceAt] .== a,Symbol("K$k")])+sum(ϕ_links_recv[ϕ_links_recv[:sinkAt] .== a,Symbol("K$k")]))*(2*(N-1))/(in_degrees[a]+out_degrees[a])
        end
    end
    # writedlm("TestPrints/gamma_$iteration.csv", γ)
    s0 = zeros(Float64, K)
    for k in 1:K
        for m in 1:links_total
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
    # writedlm("TestPrints/tau_$iteration.csv", τ)
    Elog_Θ[:] = Elog_Dirichlet(γ)
    Elog_β[:] = Elog_Beta(τ)
    # writedlm("TestPrints/ELOG_THETA_$iteration.csv",Elog_Θ)
    # writedlm("TestPrints/ELOG_BETA_$iteration.csv",Elog_β)
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
#################################################
#################################################
#################################################


end
