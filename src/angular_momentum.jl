export test_clebsch_ls

function factorial_half(n)
    @assert n ≥ 0
    s = 1.0
    for i in 2:div(n, 2)
        s *= i
    end
    return s
end

"""
    clebsch(j₁, m₁, j₂, m₂, j, m)
Calculate Clebsch-Gordan coefficients.
All spins are expressed as double their actual values.
"""
function clebsch(j₁, m₁, j₂, m₂, j, m)
    @assert abs(m₁) ≤ j₁ && iseven(j₁ - m₁)
    @assert abs(m₂) ≤ j₂ && iseven(j₂ - m₂)
    @assert abs(m) ≤ j && iseven(j - m)
    
    # j = |j₁ - j₂|, |j₁ - j₂|+1, ..., j₁ + j₂
    if !(abs(j₁-j₂) ≤ j ≤ j₁+j₂) || isodd(j-j₁-j₂) || m ≠ m₁+m₂
        return 0.0
    end
    
    # f(n) = (n/2)!
    f = factorial_half
    
    con = (j+1)*f(j₁+j₂-j)*f(j₁-j₂+j)*f(-j₁+j₂+j) / f(j₁+j₂+j+2)
    con *= f(j₁+m₁)*f(j₁-m₁)*f(j₂+m₂)*f(j₂-m₂)*f(j+m)*f(j-m)
    con = sqrt(con)
    
    sum = 0.0
    for z in max(0, -j+j₂-m₁, -j+j₁+m₂): 2: min(j₁+j₂-j, j₁-m₁, j₂+m₂)
        sum += (-1)^div(z,2)/( f(z)*f(j₁+j₂-j-z)*f(j₁-m₁-z)*f(j₂+m₂-z)*f(j-j₂+m₁+z)*f(j-j₁-m₂+z) )
    end
    
    return con*sum
end




"""
    wigner3j(j₁, j₂, j₃, m₁, m₂, m₃)
Calculate wigner 3j symbol.
All spins are expressed as double their actual values.
"""
function wigner3j(j₁, j₂, j₃, m₁, m₂, m₃)
    return (-1)^div(j₁-j₂-m₃,2)/sqrt(j₃+1)*clebsch(j₁, m₁, j₂, m₂, j₃, -m₃)
end


"""
    calc_angular_matrix_element(l₁, j₁, m₁, λ, μ, l₂, j₂, m₂)

l₁, l₂, λ, μ: actual values
j₁, m₁, j₂, m₂: double their actual values
"""
function calc_angular_matrix_element(l₁, j₁, m₁, λ, μ, l₂, j₂, m₂)
    if isodd(l₁+l₂+λ)
        return 0.0
    end
    
    M = clebsch(j₁, m₁, j₂, -m₂, 2λ, 2μ)/sqrt(2λ+1) *
    sqrt((j₁+1)*(j₂+1)/4π) * clebsch(j₁, 1, j₂, -1, 2λ, 0) 

    if isodd(div(j₂-m₂,2)) 
        M *= -1 
    end

    if isodd(div(j₁-1,2)+λ)
        M *= -1 
    end
    
    return M
end


function calc_reduced_angular_matrix_element(l₁, j₁, λ, l₂, j₂)
    if isodd(l₁+l₂+λ)
        return 0.0
    end
    
    M = sqrt((j₁+1)*(j₂+1)/4π) * clebsch(j₁, 1, j₂, -1, 2λ, 0) 

    if isodd(div(j₁-1,2)+λ)
        M *= -1 
    end
    
    return M
end




"""
    clebsch_ls(l, j, m, ms)
Calculate (l, m-ms, 1/2, ms | j, m).
All spins other than l are expressed as double their actual values.
"""
function clebsch_ls(l, j, m, ms)
    @assert j === 2l + 1 || j === 2l - 1
    @assert abs(m) ≤ j && iseven(j - m)
    @assert ms === 1 || ms === -1
    @assert abs(m-ms) ≤ 2l

    if j === 2l + 1
        if ms === +1
            return sqrt((j+m)/2j)
        end
        if ms === -1
            return sqrt((j-m)/2j)
        end
    end

    if j === 2l - 1 
        if ms === +1
            return -sqrt((j+2-m)/2(j+2))
        end
        if ms === -1
            return +sqrt((j+2+m)/2(j+2))
        end
    end
end


function test_clebsch_ls(lmax)
    for l in 0:lmax, j in 2l+1: -2: max(2l-1, 0), m in j: -2: -j
        for ms in 1: -2: -1
            if abs(m-ms) > j
                continue 
            end 
            error = clebsch(2l, m-ms, 1, ms, j, m) - clebsch_ls(l, j, m, ms)
            #println("")
            #@show l, j, m, ms
            #@show error
        end
    end

    @time for l in 0:lmax, j in 2l+1: -2: max(2l-1, 0), m in j: -2: -j
        for ms in 1: -2: -1
            if abs(m-ms) > j
                continue 
            end 
            clebsch(2l, m-ms, 1, ms, j, m)
        end
    end

    @time for l in 0:lmax, j in 2l+1: -2: max(2l-1, 0), m in j: -2: -j
        for ms in 1: -2: -1
            if abs(m-ms) > j
                continue 
            end 
            clebsch_ls(l, j, m, ms)
        end
    end
end

