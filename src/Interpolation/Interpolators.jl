module Interpolators
    abstract type Interpolator end

    function idx_function_linear_grid(x::Number, x_arr)::Int
        spacing::Float64 = (x_arr[end] - x_arr[1]) / (length(x_arr) - 1)
        return ((x - x_arr[1]) ÷ spacing) + 1   # This is an integer division
    end

    function idx_function_general(x::Number, x_arr)::Int
        for (i, val) in enumerate(x_arr)
            if (val > x)
                return i - 1
            end
        end

        return nothing
    end

    function idx_function_log_grid(x::Number, x_arr)::Int
        N = length(x_arr)
        i = floor((N - 1) * (log(x) - log(x_arr[1])) / (log(x_arr[end]) - log(x_arr[1]))) + 1
        return i
    end

    function idx_function_centered_logspaced_grid(x::Number, x_arr)
        # Use logarithmic spacing functions for each of the 2 subarrays
        center_idx = div(length(x_arr), 2) + 1

        center_element = x_arr[center_idx]
        if x > center_element
            i = idx_function_log_grid(x, x_arr[center_idx:end]) + center_idx - 1
            return i
        elseif x < center_element
            i = idx_function_log_grid(x, x_arr[1:center_idx])
            return i
        else
            return center_idx
        end
    end

    include("LinearInterpolation.jl")
    include("CSplineInterpolation.jl")
    include("ChebyInterpolation.jl")
    include("ComplexCSplineLinInterpolator.jl")
end