export legendre, test_legendre

"""
legendre(l, m, x)

Calculate legendre function
See http://www.natural-science.or.jp/article/20111227231948.php
"""
function legendre(l, m, x)
    if m < 0
        return (-1)^m * factorial(l+m)/factorial(l-m) * legendre(l, -m, x)
    end
    
    P₀ = zero(x) # P(m-1, m)
    if l < m
        return P₀
    end
    
    P₁ = sqrt((1 - x^2)^m) # P(m, m)
    for i in m+1:2m
        P₁ *= i/2
    end
    if l == m
        return P₁
    end
    
    P₂ = zero(x) 
    for i in m+1:l
        P₂ = (2i - 1)/(i - m)*x*P₁ - (i + m - 1)/(i - m)*P₀
        P₀ = P₁
        P₁ = P₂
    end
    return P₂
end


function test_legendre(l, m, legendre_exact)
    f(x) = legendre(l, m, x)
    
    p = plot(xlim=(-1,1))
    plot!(p, f)
    plot!(p, legendre_exact)
    display(p)

    xs = range(-1, 1, step=0.1)
    fs = f.(xs)
    gs = legendre_exact.(xs)
    error = sum(@. (fs-gs)^2)
    @show error
    return 
end



function spherical_harmonics(l, m, θ, φ)
    sqrt((2l+1)/4π * factorial(l-m)/factorial(l+m)) * 
    legendre(l, m, cos(θ)) * exp(im*m*φ)
end