
export test_calc_BE1_strength

function calc_dipole_matrix_element(param, spstates, n₁, n₂)
    @unpack Z, A, Nr, rs, Δr = param 
    @unpack nstates, ψs, spEs, qnums = spstates 

    l₁ = qnums[n₁].l 
    j₁ = qnums[n₁].j 

    l₂ = qnums[n₂].l 
    j₂ = qnums[n₂].j 

    if iseven(l₁+l₂) 
        return 0.0 
    end

    ME_rad = 0.0 
    for ir in 1:Nr 
        ME_rad += ψs[ir, n₁] * rs[ir] * ψs[ir, n₂]  
    end
    ME_rad *= Δr 

    ME_ang = (-1)^div(j₁+1,2) * clebsch(j₂, 1, j₁, -1, 2, 0) * sqrt((j₁+1)*(j₂+1)/4π)
    #ME_ang = calc_reduced_angular_matrix_element(l₁, j₁, 1, l₂, j₂) 

    return -(Z/A) * ME_rad * ME_ang 
end



function calc_BE1_strength(param, spstates, coeff_gs, coeff_excited)  
    @unpack Z, A, Nr, rs, Δr, Emax, lmax = param 
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    J_gs = 0 # actual value 
    M_gs = 0 # actual value 
    dim_gs = length(coeff_gs) 

    J_excited = 1 # actual value 
    M_excited = 0 # actual value 
    dim_excited = length(coeff_excited) 

    ME = 0.0 

    n₃₄ = 0 # ground states 
    for n₄ in 1:nstates
        if occ[n₄] == 1.0
            continue 
        end 
        spE₄ = spEs[n₄] 
        l₄ = qnums[n₄].l 
        j₄ = qnums[n₄].j 

        for n₃ in 1:n₄ # n₃ ≤ n₄
            if occ[n₃] == 1.0
                continue 
            end
            spE₃ = spEs[n₃]
            l₃ = qnums[n₃].l 
            j₃ = qnums[n₃].j 

            if isodd(l₃+l₄+J_gs) || !(abs(j₃-j₄) ≤ 2J_gs ≤ j₃+j₄) || isodd(2J_gs-j₃-j₄)
                continue 
            end
            if spE₃ + spE₄ > Emax 
                continue 
            end 

            n₃₄ += 1 

            n₁₂ = 0 # excited states
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
        
                    if isodd(l₁+l₂+J_excited) || !(abs(j₁-j₂) ≤ 2J_excited ≤ j₁+j₂) || isodd(2J_excited-j₁-j₂)
                        continue 
                    end 
                    if spE₁ + spE₂ > Emax 
                        continue 
                    end
        
                    n₁₂ += 1

                    temp = coeff_excited[n₁₂] * coeff_gs[n₃₄] * 
                    1/sqrt(3(j₃+1) * (1+(n₁===n₂)) * (1+(n₃===n₄)))

                    if n₂ === n₄ 
                        ME += temp * (-1)^div(j₁+j₃,2)
                        calc_dipole_matrix_element(param, spstates, n₁, n₃) 
                    end

                    if n₁ === n₃ 
                        ME += temp * #(-1)^div(j₂+j₄,2) *
                        calc_dipole_matrix_element(param, spstates, n₂, n₄) 
                    end

                    if n₂ === n₃ 
                        ME += temp * (-1)^div(j₁+j₄,2)
                        calc_dipole_matrix_element(param, spstates, n₁, n₄) 
                    end 

                    if n₁ === n₄ 
                        ME += temp * #(-1)^div(j₂+j₃,2) * 
                        calc_dipole_matrix_element(param, spstates, n₂, n₃) 
                    end

                end
            end
            @assert n₁₂ === dim_excited
        end
    end
    @assert n₃₄ === dim_gs 

    BE1 = ME^2 

    return BE1 
end

function test_calc_BE1_strength(param; howmany=50, Γ=0.2, figname="test")
    @unpack Emax, lmax = param 

    spstates = calc_single_particle_states(param)
    calc_occ!(spstates, param)

    # ground states 
    J_gs = 0 
    M_gs = 0 
    Hmat_3body = make_three_body_Hamiltonian(param, spstates, J_gs) 
    @time Es_gs, coeffs_gs = eigen(Hmat_3body)
    E_gs = Es_gs[1]
    coeff_gs = coeffs_gs[:, 1] 
    @show E_gs 
    
    J_excited = 1 
    Hmat_3body = make_three_body_Hamiltonian(param, spstates, J_excited) 
    @time Es_excited, coeffs_excited = eigen(Hmat_3body)
    @show Es_excited[1:2]

    BE1s = zeros(Float64, length(Es_excited))
    @time for k in 1:length(Es_excited)
        BE1s[k] = calc_BE1_strength(param, spstates, coeff_gs, coeffs_excited[:,k]) 
    end
    @show BE1s[1:2]

    Es = range(0, 5.0,  step=0.01) 
    fs = zeros(Float64, length(Es)) 
    for k in 1:length(Es_excited) 
        for iE in 1:length(Es)
            E = Es[iE] 
            fs[iE] += 3 * (Γ/π) * 1/((E - Es_excited[k] + E_gs)^2 + Γ^2) * BE1s[k]
        end
    end
    p = plot(xlabel="E [MeV]", ylabel="B(E1)", title="Emax=$(Emax)[MeV]  lmax=$(lmax)  Γ=$Γ[MeV]", ylim=(0,5))
    plot!(p, Es, fs; label="total")
    #savefig("./Figure/" * figname * ".png")
    display(p)
end