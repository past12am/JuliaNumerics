module LinearInterpolation
    
    # If I want to interpolate a complex function --> 2D interpolation/separately for real/imag part?
    struct LinearInterpolator 
        x::AbstractArray{<:Number}    # TODO This should be a pointer?
        f::AbstractArray{<:Number}    # TODO This should be a pointer?
    end

    # Spacing function yields how the grid is interpolated, it should have signature
    #   idx_function(x, num_points, lower_bound, upper_bound) and return an index i s.t. (x_i < x < x_i+1)
    function interpolate(x::Number, lin_interp::LinearInterpolator, idx_function::Function)
        if (lin_interp.x[1] < x && x < lin_interp.x[end])
            i::Int = idx_function(x, lin_interp.x)
            return lin_interp.f[i] + (lin_interp.f[i+1] - lin_interp.f[i]) / (lin_interp.x[i+1] - lin_interp.x[i]) * (x - lin_interp.x[i])
        
        elseif (x == lin_interp.x[1])
            return lin_interp.f[1]

        elseif (x == lin_interp.x[end])
            return lin_interp.f[end]
        end

        throw(ErrorException("Out of bounds for interpolation"))
    end
end