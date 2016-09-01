module genr_network
using Distributions
using Devectorize
N=75
K=5
α=Array(Float64, K)
for k in 1:K
    α[k] = 0.01
end
ϵ = 1e-30
dbeta=Beta()
β_diag=
end
