module SphericalWoodsSaxon

using Plots
using LinearAlgebra
using Parameters
using KrylovKit
using ProgressMeter

include("./single_particle_states.jl")
include("./angular_momentum.jl")
include("./legendre.jl")

include("./three_body_Hamiltonian.jl")
include("./three_body_density.jl")
include("./three_body.jl")

export PhysicalParam

@with_kw struct PhysicalParam{T} @deftype Float64
    ħc = 197.
    mc² = 938.
    M = ħc^2/2mc²
    
    # particle number
    Z::Int64 = 6; @assert iseven(Z) 
    N::Int64 = 14; @assert iseven(N) 
    A::Int64 = Z + N; @assert A === Z + N
    
    # parameters of Woods-Saxon potential
    V₀ = -38.76*0.963 # [MeV]
    V₁ = -25.63*0.963
    r₀ = 1.25 # [fm]
    R₀ = r₀*A^(1/3) # [fm]
    a = 0.65 # [fm]
    #κ = 0.44
    
    # model space
    Emax = 5 # [MeV]
    lmax::Int64 = 5
    Λmax::Int64 = 2lmax+1; @assert isodd(Λmax)

    # parameters of neutron-neutron interaction 
    a_nn = -15.0 # [fm]
    v₀_nn = 2π^2*ħc^2/mc² * 2a_nn/(π - 2a_nn*sqrt(mc²*Emax/ħc^2))
    v_rho = -v₀_nn 
    R_rho = r₀*A^(1/3)
    a_rho = 0.67
    
    # radial mesh
    Nr::Int64 = 100
    Δr = 0.2
    rs::T = range(Δr, Nr*Δr, length=Nr)
end

end # module
