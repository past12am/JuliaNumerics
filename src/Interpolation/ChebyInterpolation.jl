module ChebyInterpolation
    import ....Numerics.ChebyshevPolynomials as ChPoly
    
    struct ChebyInterpolatorFirstKind
        #linear_interpolator::LinearInterpolator
        d::AbstractArray{<:Number}
        M::Int

        ChebyInterpolatorFirstKind(x::AbstractArray, f::AbstractArray, M::Int) = begin
            return new(calc_di_first_kind(x, f, M, length(x)), M)
        end
    end

    function calc_di_first_kind(x::AbstractArray, f::AbstractArray, M::Int, N::Int)::AbstractArray
        d = zeros(Float64, M)

        for i = 1:M
            for k = 1:N
                d[i] += f[k] * ChPoly.chebyshev_poly_first_kind(x[k], i-1)     # 1-based indexing is a mess
            end 

            d[i] *= 2.0 / N
        end

        return d
    end

    function interpolate(x::Number, interpolator::ChebyInterpolatorFirstKind)
        #res = interpolator.d[1] / 2.0
        #for i = 2:interpolator.M
        #    res += interpolator.d[i] * ChPoly.chebyshev_poly_first_kind(x, i-1)
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



    struct ChebyInterpolatorSecondKind
        #linear_interpolator::LinearInterpolator
        c::AbstractArray{<:Number}
        n::Int

        fType::Type

        ChebyInterpolatorSecondKind(x::AbstractArray, f::AbstractArray) = begin
            n = length(x)
            return new(calc_ci_second_kind(x, f, n), n, eltype(f))
        end
    end

    function calc_ci_second_kind(x::AbstractArray, f::AbstractArray, n::Int)::AbstractArray
        c = zeros(eltype(f), n)

        for i = 0:n-1
            for j = 0:n-1
                c[i+1] += 2.0 / ((n - 1) + 2.0) * (1.0 - x[j+1]^2) * f[j+1] * ChPoly.chebyshev_poly_second_kind_cosTheta(x[j+1], i)
            end
        end

        return c
    end

    function interpolate(x::Number, interpolator::ChebyInterpolatorSecondKind)
        res = zero(interpolator.fType)

        for i = 0:interpolator.n-1
            res += interpolator.c[i+1] * ChPoly.chebyshev_poly_second_kind_cosTheta(x, i)
        end

        return res
    end
end