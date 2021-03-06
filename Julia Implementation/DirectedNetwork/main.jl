module Main

include("../DirectedNetwork/network_genr.jl")
include("../DirectedNetwork/youtilities.jl")
include("../DirectedNetwork/inference.jl")

using Plots

plot(1:length(Inference.ELBO),Inference.ELBO)
for k in 1:NetworkGenr.K
    println(Inference.τ[1,k]/(sum(Inference.τ[:,k])))
end

[NetworkGenr.β[i,i] for i in 1:size(NetworkGenr.β)[1]]
