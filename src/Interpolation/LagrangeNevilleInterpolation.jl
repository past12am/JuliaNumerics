module LagrangeNevilleInterpolation

    function interpolate(x_vals::AbstractVector{<:Number}, y_vals::AbstractVector{<:Number}, x::Number)
        M = length(y_vals)
        T = float(promote_type(eltype(x_vals), eltype(y_vals), typeof(x)))

        C = zeros(T, M)
        D = zeros(T, M)

        return interpolate!(C, D, x_vals, y_vals, x)
    end

    function interpolate!(C::AbstractVector{T}, D::AbstractVector{T}, x_vals::AbstractVector{<:Number}, y_vals::AbstractVector{<:Number}, x) where {T<:Number}
        
        M = length(y_vals)

        # Init
        #   Initial Approximation
        closest_idx::Int = 1
        min_diff = abs(x_vals[1] - x)

        #   Init grids
        for i = 1:M
            C[i] = y_vals[i]
            D[i] = y_vals[i]

            diff = abs(x_vals[i] - x)
            if (diff < min_diff)
                closest_idx = i
                min_diff = diff
            end
        end

        #   Init y estimate
        y = C[closest_idx]
        closest_idx -= 1

        #   Recurrence
        for m = 1:M-1
            for i = 1:(M-m)
                ho = x_vals[i] - x
                hp = x_vals[i+m] - x

                den = ho - hp

                w = C[i+1] - D[i]

                if (den == 0.0)
                    throw(ArgumentError("identical nodes"))      # This error can occur only if two input xa’s are (to within roundoff) identical.
                end
                
                wden = w/den

                C[i] = ho * wden
                D[i] = hp * wden
            end

            if (2 * closest_idx < (M - m))
                y += C[closest_idx + 1]
            else
                y += D[closest_idx]
                closest_idx -= 1
            end
        end

        return y

    end

end