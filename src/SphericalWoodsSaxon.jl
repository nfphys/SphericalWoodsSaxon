module SphericalWoodsSaxon

using Plots
using LinearAlgebra
using Parameters
using KrylovKit
using ProgressMeter
using Roots

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
    
    # particle number of core nucleus
    Z::Int64 = 6; @assert iseven(Z) 
    N::Int64 = 14; @assert iseven(N) 
    A::Int64 = Z + N; @assert A === Z + N

    M = ħc^2/2mc²*(1 + 1/A)
    
    # parameters of Woods-Saxon potential
    #= O-24
    V₀ = -43.2 # [MeV]
    V₁ = 0.44*0.73*V₀
    r₀ = 1.25 # [fm]
    R₀ = r₀*A^(1/3) # [fm]
    a = 0.67 # [fm]
    V_gaus = 0.0
    μ_gaus = 0.09
    =#

    # C-22
    V₀ = -38.76 # [MeV]
    V₁ = -25.63/1.25^2
    r₀ = 1.25 # [fm]
    R₀ = r₀*A^(1/3) # [fm]
    a = 0.65 # [fm]
    V_gaus = 4.66
    μ_gaus = 0.09
    
    # model space
    Emax = 5 # [MeV]
    lmax::Int64 = 5

    # radial mesh
    Nr::Int64 = 40
    Δr = 0.5
    rs::T = range(Δr, Nr*Δr, length=Nr)
    ir_matching::Int64 = floor(Int, R₀/Δr)

    # parameters of neutron-neutron interaction
    #= O-24 
    a_nn = -15.0 # [fm]
    v₀_nn = 2π^2*ħc^2/mc² * 2a_nn/(π - 2a_nn*sqrt(mc²*Emax/ħc^2))
    v_rho = 814.2
    R_rho = R₀
    a_rho = 0.67
    =#

    # C-22
    a_nn = -15.0 # [fm]
    v₀_nn = 2π^2*ħc^2/mc² * 2a_nn/(π - 2a_nn*sqrt(mc²*Emax/ħc^2))
    v_rho = -v₀_nn 
    R_rho = R₀
    a_rho = 0.67
end 

end # module
