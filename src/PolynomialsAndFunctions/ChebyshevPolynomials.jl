module ChebyshevPolynomials    
    
    function chebyshev_poly(x::Real, N::Integer)
        if N==0
            return 1.0

        elseif N==1
            return x

        elseif N >= 2
            T_pp::Float64 = 1.0
            T_p::Float64 = x

            T_cur::Float64 = 0.0
            for i = 2:N
                T_cur = 2.0 * x * T_p - T_pp

                T_pp = T_p
                T_p = T_cur
            end

            return T_cur
        end

        throw(ErrorException("Polynomial of degree $(N) seems weird"))
    end

end