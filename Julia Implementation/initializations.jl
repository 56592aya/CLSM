module Initializations
#global


#first round init
using GenericFuncs: Elog_Dirichlet, Elog_Beta

Elog_Θ = Elog_Dirichlet(γ)
Elog_β = Elog_Beta(τ0, τ1)

export Elog_Θ, Elog_β

end
