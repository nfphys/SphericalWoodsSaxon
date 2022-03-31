export test_calc_single_particle_states

@with_kw struct QuantumNumbers @deftype Int64
    l = 0
    j = 1
end

@with_kw struct SingleParticleStates
    nstates::Int64
    ψs::Matrix{Float64}; @assert size(ψs, 2) === nstates
    spEs::Vector{Float64}; @assert length(spEs) === nstates
    qnums::Vector{QuantumNumbers}; @assert length(qnums) === nstates
    occ::Vector{Float64}; @assert length(occ) === nstates
end




function make_single_particle_Hamiltonian(param, qnum)
    @unpack M, Nr, Δr, rs, V₀, V₁, r₀, R₀, a, V_gaus, μ_gaus = param
    @unpack l, j = qnum
    
    # potential
    Vs = zeros(Float64, Nr)
    
    # central part
    @. Vs += V₀/(1 + exp((rs-R₀)/a))
    
    # spin-orbit part
    ls = (j*(j+2) - 4l*(l+1) - 3)/8
    @. Vs += V₁*ls*(r₀*r₀/(a*rs))*exp((rs-R₀)/a)/(1 + exp((rs-R₀)/a))^2

    # gaussian potential for s-wave
    if l === 0
        @. Vs += V_gaus*exp(-μ_gaus*rs*rs)
    end
    
    # centrifugal part
    @. Vs += M*l*(l+1)/rs^2
    
    dv = zeros(Float64, Nr)
    ev = zeros(Float64, Nr-1)
    
    @. dv = 2M/Δr^2 + Vs
    @. ev = -M/Δr^2
    
    return SymTridiagonal(dv, ev)
end

function calc_single_particle_states(param)
    @unpack Nr, Δr, Emax, lmax = param 
    
    ψs = zeros(Float64, Nr, Nr*(lmax+1))
    spEs = zeros(Float64, Nr*(lmax+1))
    qnums = Vector{QuantumNumbers}(undef, Nr*(lmax+1))
    occ = zeros(Float64, Nr*(lmax+1))
    
    nstates = 0
    for l in 0:lmax, j in (2l+1): -2: max(2l-1, 0)
        qnum = QuantumNumbers(l=l,j=j)
        
        Hmat = make_single_particle_Hamiltonian(param, qnum)
        vals, vecs = eigen(Hmat)
        #@show dot(vecs[:,1], vecs[:,1])
        
        @. vecs /= sqrt(Δr)
        
        #@show dot(vecs[:,1], vecs[:,1])*Δr
        
        for i in 1:Nr
            if vals[i] > Emax
                break
            end
            nstates += 1
            ψs[:,nstates] = vecs[:,i]
            spEs[nstates] = vals[i]
            qnums[nstates] = qnum
        end
    end
    
    p = sortperm(spEs[1:nstates])
    
    SingleParticleStates(nstates, ψs[:,p], spEs[p], qnums[p], occ[p])
end


function calc_occ!(spstates, param)
    @unpack N = param
    @unpack nstates, qnums, occ = spstates
    
    fill!(occ, 0)
    n = 0
    for i in 1:nstates
        @unpack l, j = qnums[i]
        if n + (j+1) ≤ N
            occ[i] = 1
            n += (j+1)
        elseif n < N
            occ[i] = (N - n)/(j+1)
            n = N
        end
    end
    @assert n === N
    return
end

function show_spstates(spstates; Emax_show=0.0)
    @unpack nstates, spEs, qnums, occ = spstates 
    println("")
    for i in 1:nstates
        if spEs[i] > Emax_show
            break 
        end
        println("i = ", i, ": ")
        @show spEs[i] occ[i] qnums[i]
    end
end


function test_calc_single_particle_states(param; Emax_show=2.0)
    @time spstates = calc_single_particle_states(param)
    calc_occ!(spstates, param)
    show_spstates(spstates; Emax_show=Emax_show)
end