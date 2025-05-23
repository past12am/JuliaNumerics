module LegendrePolynomials

    # Recursive Definition (do not use, slow)
    function legendre_poly_rec(x::Number, N::Integer)
        if N >= 2
            return ((2.0 * N - 1.0) * x * legendre_poly(x, N-1) - (N - 1.0) * legendre_poly(x, N-2)) / N

        elseif N == 1
            return x

        elseif N == 0
            return 1.0
        end

        throw(ErrorException("N=$(N) is an odd Legendre Polynomial request"))
    end

    function derivative_legendre_poly_rec(x::Number, N::Integer)
        if(N > 1)
            return N * (legendre_poly(x, N-1) - x * legendre_poly(x, N)) / (1.0 - x^2)
        end
        
        return 0
    end

    

    # Looped Definition
    function legendre_poly_and_deriv(x::Real, N::Integer)
        if N == 0
            return (1.0, 0.0)
        elseif N == 1
            return (x, 1.0)
        end
    
        p_p = 1.0       # p_previous
        p = x           # p
        dp = 1.0        # dp
    
        for i in 2:N
            # Calculate Legendre Polynomial
            p_n = ((2.0 * i - 1.0) * x * p - (i - 1.0) * p_p) / i
    
            # Calculate Derivative of Legendre Polynomial
            dp_n = N * (x * p_n - p) / (x^2 - 1.0)
    
            p_p = p
            p = p_n
            dp = dp_n
        end
    
        return (p, dp)
    end
    
    function legendre_poly(x::Real, N::Integer)
        return legendre_poly_and_deriv(x, N)[1]
    end
end