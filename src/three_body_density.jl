
export test_calc_two_particle_density, calc_two_particle_density


#=
function calc_two_particle_density(param)
    @unpack Nr, rs, Δr, Emax, lmax, R₀ = param 
    @show R₀

    spstates = calc_single_particle_states(param)
    calc_occ!(spstates, param)
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    J = 0
    M = 0
    Hmat_3body = make_three_body_Hamiltonian(param, spstates, J)

    @time vals, vecs, info = eigsolve(Hmat_3body, 1, :SR, eltype(Hmat_3body))
    E_gs = vals[1]
    coeff_gs = vecs[1]

    r = R₀
    ir = floor(Int, r/Δr)
    θs = range(0, π, length=100+1)
    Nθ = length(θs)
    ρ₂ = zeros(Float64, Nθ) 

    σ₁ = +1
    σ₂ = -1

    n₁₂ = 0
    for n₂ in 1:nstates
        if occ[n₂] == 1.0
            continue 
        end
        spE₂ = spEs[n₂]
        l₂ = qnums[n₂].l
        j₂ = qnums[n₂].j

        for n₁ in 1:n₂ # n₁ ≤ n₂
            if occ[n₁] == 1.0
                continue 
            end
            spE₁ = spEs[n₁]
            l₁ = qnums[n₁].l
            j₁ = qnums[n₁].j

            if isodd(l₁+l₂+J) || !(abs(j₁-j₂) ≤ 2J ≤ j₁+j₂) || isodd(2J-j₁-j₂)
                continue 
            end 
            if spE₁ + spE₂ > Emax 
                continue 
            end

            n₁₂ += 1

            @assert l₁ === l₂
            @assert j₁ === j₂
            
            
            @. ρ₂ += coeff_gs[n₁₂] * 
            0.5 * (-1)^l₁ * sqrt((j₁+1)/(2l₁+1)) * 
            ψs[ir, n₁]/rs[ir] * ψs[ir, n₂]/rs[ir] * 
            sqrt((2l₁+1)/4π) * legendre(l₁, 0, cos(θs)) 

            #=
            @. ρ₂ += coeff_gs[n₁₂] * 1/sqrt(1 + (n₁===n₂)) *
            clebsch(j₁, 2M-σ₂, j₂, σ₂, J, M) * 
            clebsch_ls(l₁, j₁, 2M-σ₂, σ₁) * 
            =#
            
            #=
            @. ρ₂ += ψs[ir, n₁]/rs[ir] * ψs[ir, n₂]/rs[ir] *
            clebsch(j₁, 2M-σ₂, j₂, σ₂, 2J, 2M) * 
            clebsch_ls(l₁, j₁, 2M-σ₂, σ₁) * Y(l₁, div(2M-σ₁-σ₂,2), θ, φ) *
            clebsch_ls(l₂, j₂,    σ₂, σ₂) * sqrt((2l₂+1)/4π) 
            =#

        end
    end

    @. ρ₂ = ρ₂^2 /4π

    p = plot(xlabel="θ/π", ylabel="two-particle density [fm⁻⁶]", title="Spherical")
    plot!(p, θs/π, ρ₂)
    display(p)
end
=#



function calc_uncorrelated_2pwf(param, spstates, J, M, n₁, n₂, ir, θ, φ, σ₁, σ₂)
    @unpack Nr, rs, Δr, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    #ir = floor(Int, r/Δr)
    @assert 1 ≤ ir ≤ Nr
    r = rs[ir]

    Y = spherical_harmonics

    l₁ = qnums[n₁].l
    j₁ = qnums[n₁].j

    l₂ = qnums[n₂].l
    j₂ = qnums[n₂].j

    if abs(2M-σ₂) > j₁ || abs(2M-σ₁-σ₂) > 2l₁
        return 0.0 
    end

    
    a = 0.5 * (-1)^l₁ * sqrt((j₁+1)/(2l₁+1))
    b = sqrt(2l₁+1) * clebsch(j₁, 2M-σ₂, j₂, σ₂, 2J, 2M) * 
    clebsch_ls(l₁, j₁, 2M-σ₂, σ₁) * clebsch_ls(l₂, j₂, σ₂, σ₂) 
    @assert a ≈ b
    

    ψ₂ = ψs[ir, n₁]/rs[ir] * ψs[ir, n₂]/rs[ir] *
    clebsch(j₁, 2M-σ₂, j₂, σ₂, 2J, 2M) * 
    clebsch_ls(l₁, j₁, 2M-σ₂, σ₁) * Y(l₁, div(2M-σ₁-σ₂,2), θ, φ) *
    clebsch_ls(l₂, j₂,    σ₂, σ₂) * sqrt((2l₂+1)/4π) 
end

function calc_two_particle_density(param, spstates, J, M, coeff, ir, θ, φ, σ₁, σ₂)
    @unpack Nr, rs, Δr, Emax, lmax, R₀ = param 
    @unpack nstates, ψs, spEs, qnums, occ = spstates 
    #@show R₀

    ψ₂ = 0.0 + 0.0im 

    n₁₂ = 0
    for n₂ in 1:nstates
        if occ[n₂] == 1.0
            continue 
        end
        spE₂ = spEs[n₂]
        l₂ = qnums[n₂].l
        j₂ = qnums[n₂].j

        for n₁ in 1:n₂ # n₁ ≤ n₂
            if occ[n₁] == 1.0
                continue 
            end
            spE₁ = spEs[n₁]
            l₁ = qnums[n₁].l
            j₁ = qnums[n₁].j

            if isodd(l₁+l₂+J) || !(abs(j₁-j₂) ≤ 2J ≤ j₁+j₂) || isodd(2J-j₁-j₂)
                continue 
            end 
            if spE₁ + spE₂ > Emax 
                continue 
            end

            n₁₂ += 1

            if J === 0
                @assert l₁ === l₂
                @assert j₁ === j₂
            end
            
            temp = 0.0
            temp += calc_uncorrelated_2pwf(param, spstates, J, M, n₁, n₂, ir, θ, φ, σ₁, σ₂)/√2
            temp += calc_uncorrelated_2pwf(param, spstates, J, M, n₂, n₁, ir, θ, φ, σ₁, σ₂)/√2 * (-1)^(j₁-j₂-J)
            if n₁ === n₂ 
                temp /= √2
            end

            ψ₂ += coeff[n₁₂] * temp

        end
    end
    @assert n₁₂ === length(coeff)

    ρ₂ = abs2(ψ₂)
end

function test_calc_two_particle_density(param; J=0, M=0, β=0.0, σ₁=1, σ₂=-1, φ=0.0)
    @unpack Nr, rs, Δr, Emax, lmax, R₀ = param 
    @show R₀

    spstates = calc_single_particle_states(param)
    calc_occ!(spstates, param)
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    Hmat_3body = make_three_body_Hamiltonian(param, spstates, J)

    @time vals, vecs, info = eigsolve(Hmat_3body, 1, :SR, eltype(Hmat_3body))
    coeff = vecs[1]

    #r = R₀
    θs = range(0, π, length=100+1)
    Nθ = length(θs)
    φ = 0

    ρ₂ = zeros(Float64, Nr, Nθ)
    @time for iθ in 1:Nθ, ir in 1:Nr 
        r = rs[ir]
        θ = θs[iθ]
        ρ₂[ir, iθ] = 4π*r^2 * 4π*r^2 * sin(θ) * 
        calc_two_particle_density(param, spstates, J, M, coeff, ir, θ, φ, σ₁, σ₂)
    end

    #=
    p = plot(xlabel="r [fm]", ylabel="θ/π")
    heatmap!(p, rs, θs/π, ρ₂')
    display(p)
    =#

    f(θ) = calc_two_particle_density(param, spstates, J, M, coeff, ir, θ, φ, σ₁, σ₂)
    @time ρ₂ = f.(θs)

    p = plot(xlabel="θ/π", ylabel="two-particle density [fm⁻⁶]", title="Spherical")
    plot!(p, θs/π, ρ₂)
    display(p)
end





#=
function calc_single_particle_density(param)
    @unpack Nr, rs, Δr, Emax, lmax, R₀ = param 
    @show R₀

    spstates = calc_single_particle_states(param)
    calc_occ!(spstates, param)
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    J = 0
    Hmat_3body = make_three_body_Hamiltonian(param, spstates, J)

    @time vals, vecs, info = eigsolve(Hmat_3body, 1, :SR, eltype(Hmat_3body))
    E_gs = vals[1]
    coeff_gs = vecs[1]

    # calculate single-particle density
    ρ₁ = zeros(Float64, Nr) 

    n₃₄ = 0
    @time for n₄ in 1:nstates
        if occ[n₄] == 1.0
            continue 
        end 
        spE₄ = spEs[n₄] 
        l₄ = qnums[n₄].l 
        j₄ = qnums[n₄].j 

        for n₃ in 1:nstates # n₃ ≤ n₄
            if occ[n₃] == 1.0
                continue 
            end
            spE₃ = spEs[n₃]
            l₃ = qnums[n₃].l 
            j₃ = qnums[n₃].j 

            if isodd(l₃+l₄+J) || !(abs(j₃-j₄) ≤ 2J ≤ j₃+j₄) || isodd(2J-j₃-j₄)
                continue 
            end
            if spE₃ + spE₄ > Emax 
                continue 
            end 

            n₃₄ += 1 

            n₁₂ = 0
            for n₂ in 1:nstates
                if occ[n₂] == 1.0
                    continue 
                end
                spE₂ = spEs[n₂]
                l₂ = qnums[n₂].l
                j₂ = qnums[n₂].j
        
                for n₁ in 1:nstates # n₁ ≤ n₂
                    if occ[n₁] == 1.0
                        continue 
                    end
                    spE₁ = spEs[n₁]
                    l₁ = qnums[n₁].l
                    j₁ = qnums[n₁].j
        
                    if isodd(l₁+l₂+J) || !(abs(j₁-j₂) ≤ 2J ≤ j₁+j₂) || isodd(2J-j₁-j₂)
                        continue 
                    end 
                    if spE₁ + spE₂ > Emax 
                        continue 
                    end
        
                    n₁₂ += 1
                    
                    if n₂ === n₄ 
                        @views @. ρ₁ += (1/4π)  * coeff_gs[n₁₂] * coeff_gs[n₃₄] *
                        ψs[:, n₁] * ψs[:, n₃] / rs^2 
                    end

                end
            end
        end
    end

    num_particle = sum(@. rs^2*ρ₁)*4π*Δr
    @show num_particle

    p = plot(xlabel="r [fm]", ylabel="single-particle density [fm⁻³]")
    plot!(p, rs, ρ₁)
    display(p)
end
=#