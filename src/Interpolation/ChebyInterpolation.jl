module ChebyInterpolation
    import ....Numerics.ChebyshevPolynomials as ChPoly
    
    struct ChebyInterpolator
        #linear_interpolator::LinearInterpolator
        d::AbstractArray{<:Number}
        M::Int

        ChebyInterpolator(x::AbstractArray, f::AbstractArray, M::Int) = begin
            return new(calc_di(x, f, M, length(x)), M)
        end
    end

    function calc_di(x::AbstractArray, f::AbstractArray, M::Int, N::Int)::AbstractArray
        d = zeros(Float64, M)

        for i = 1:M
            for k = 1:N
                d[i] += f[k] * ChPoly.chebyshev_poly(x[k], i-1)     # 1-based indexing is a mess
            end 

            d[i] *= 2.0 / N
        end

        return d
    end

    function interpolate(x::Number, interpolator::ChebyInterpolator)
        #res = interpolator.d[1] / 2.0
        #for i = 2:interpolator.M
        #    res += interpolator.d[i] * ChPoly.chebyshev_poly(x, i-1)
        #end
        #return res

        e_nn::Float64 = 0.0
        e_n::Float64 = 0.0

        e::Float64 = 0.0
        
        for j = interpolator.M:-1:2
            e = 2.0 * x * e_n - e_nn + interpolator.d[j]

            e_nn = e_n
            e_n = e
        end

        return x * e_n - e_nn + interpolator.d[1] / 2.0
    end
end