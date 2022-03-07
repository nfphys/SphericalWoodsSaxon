
export test_make_three_body_Hamiltonian


function calc_Vnn_matrix_element(param, spstates, J, Vnn, n₁, n₂, n₃, n₄)
    @unpack Nr, rs, Δr, Emax, lmax, v₀_nn, v_rho, R_rho, a_rho = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    l₁ = qnums[n₁].l
    j₁ = qnums[n₁].j

    l₂ = qnums[n₂].l
    j₂ = qnums[n₂].j

    l₃ = qnums[n₃].l 
    j₃ = qnums[n₃].j 

    l₄ = qnums[n₄].l 
    j₄ = qnums[n₄].j 

    ME_rad = 0.0
    for ir in 1:Nr 
        ME_rad += ψs[ir, n₁] * ψs[ir, n₂] * Vnn[ir] *
                    ψs[ir, n₃] * ψs[ir, n₄] /rs[ir]^2
    end
    ME_rad *= Δr

    
    if isodd(l₁+l₂+J) || isodd(l₃+l₄+J)
        return 0.0 
    end

    Vnn = (-1)^(l₁+l₃) * sqrt((j₁+1)*(j₃+1))/4π * 
    clebsch(j₁, 1, 2J, 0, j₂, 1) * clebsch(j₃, 1, 2J, 0, j₄, 1) * ME_rad

    if n₁ === n₂ 
        Vnn /= √2
    end
    if n₃ === n₄
        Vnn /= √2 
    end

    return Vnn

end



function make_three_body_Hamiltonian(param, spstates, J)
    @unpack Nr, rs, Δr, Emax, lmax, v₀_nn, v_rho, R_rho, a_rho = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    # calculate the size of the three-body Hamiltonian
    dim = 0
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

            dim += 1
            #@show l₁, l₂
        end
    end
    @show Emax lmax nstates dim 

    Hmat_3body = zeros(Float64, dim, dim)

    Vnn = zeros(Float64, Nr)
    @. Vnn = v₀_nn + v_rho/(1 + exp((rs - R_rho)/a_rho))

    n₃₄ = 0
    @time for n₄ in 1:nstates
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

            if isodd(l₃+l₄+J) || !(abs(j₃-j₄) ≤ 2J ≤ j₃+j₄) || isodd(2J-j₃-j₄)
                continue 
            end
            if spE₃ + spE₄ > Emax 
                continue 
            end 

            n₃₄ += 1 

            Hmat_3body[n₃₄, n₃₄] += spE₃ + spE₄

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
                    if n₁₂ > n₃₄ # n₁₂ ≤ n₃₄
                        continue 
                    end

                    Hmat_3body[n₁₂, n₃₄] += 
                    calc_Vnn_matrix_element(param, spstates, J, Vnn, n₁, n₂, n₃, n₄)

                end
            end
            #@show n₁₂
        end
    end

    return Symmetric(Hmat_3body)
end


function test_make_three_body_Hamiltonian(param; J=0)
    spstates = calc_single_particle_states(param)
    calc_occ!(spstates, param)

    Hmat_3body = make_three_body_Hamiltonian(param, spstates, J)
    @show eigmin(Hmat_3body)
end