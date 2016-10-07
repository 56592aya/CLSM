module main
include("generic_funcs.jl")
include("genr_network.jl")
include("inference.jl")
using GenericFuncs.estimate_beta
using Inference:τ0, τ1
using GenrNetwork.β
estimate_beta(τ0, τ1)
diag(β)

end
