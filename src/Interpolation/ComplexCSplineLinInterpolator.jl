module ComplexCSplineLinInterpolation

    import ..LinearInterpolation as LinInterp
    import ..CSplineInterpolation as CSplineInterp

    struct ComplexCSplineLinInterpolator
        imaginaries::AbstractArray{<:Real}

        csplines_f_real::AbstractArray{CSplineInterp.CSplineInterpolator}
        csplines_f_imag::AbstractArray{CSplineInterp.CSplineInterpolator}

        idx_function_imaginaries::Function
        idx_function_real_grid::Function

        # We assume the imaginary part to make up the first dimension of f_real and f_imag
        ComplexCSplineLinInterpolator(x_real_arr::AbstractArray{<:Real}, x_imag_arr::AbstractArray{<:Real}, f_real_arr::AbstractArray{<:Real, 2}, f_imag_arr::AbstractArray{<:Real, 2}, idx_function_imaginaries::Function, idx_function_real_grid::Function) = begin
            
            csplines_f_real = Array{CSplineInterp.CSplineInterpolator, 1}(undef, length(x_imag_arr))
            csplines_f_imag = Array{CSplineInterp.CSplineInterpolator, 1}(undef, length(x_imag_arr))
            
            for i in eachindex(x_imag_arr)
                csplines_f_real[i] = CSplineInterp.CSplineInterpolator(x_real_arr, f_real_arr[i, :])
                csplines_f_imag[i] = CSplineInterp.CSplineInterpolator(x_real_arr, f_imag_arr[i, :])
            end

            return new(x_imag_arr, csplines_f_real, csplines_f_imag, idx_function_imaginaries, idx_function_real_grid)
        end
    end

    # Spacing function yields how the grid is interpolated, it should have signature
    #   idx_function(x, num_points, lower_bound, upper_bound) and return an index i s.t. (x_i < x < x_i+1)
    function interpolate(x::Complex, compl_spline_lin_interp::ComplexCSplineLinInterpolator)
        
        x_real = real(x)
        x_imag = imag(x)

        # Check if the requested x is within bounds TODO also for real part
        if (x_imag < compl_spline_lin_interp.imaginaries[1])
            error("requesting spline at $(x_imag) below bounds $(compl_spline_lin_interp.imaginaries[1])")
        end
        if (x_imag > compl_spline_lin_interp.imaginaries[end])
            error("requesting spline at $(x_imag) above bounds $(compl_spline_lin_interp.imaginaries[end])")
        end

        # Find idx of the spline to use
        idx = compl_spline_lin_interp.idx_function_imaginaries(x_imag, compl_spline_lin_interp.imaginaries)

        # Get real and imaginary parts at x_real < x
        f_real_prev = CSplineInterp.interpolate(x_real, compl_spline_lin_interp.csplines_f_real[idx], compl_spline_lin_interp.idx_function_real_grid)
        f_imag_prev = CSplineInterp.interpolate(x_real, compl_spline_lin_interp.csplines_f_imag[idx], compl_spline_lin_interp.idx_function_real_grid)

        # Get real and imaginary parts at x x_real
        f_real_next = CSplineInterp.interpolate(x_real, compl_spline_lin_interp.csplines_f_real[idx+1], compl_spline_lin_interp.idx_function_real_grid)
        f_imag_next = CSplineInterp.interpolate(x_real, compl_spline_lin_interp.csplines_f_imag[idx+1], compl_spline_lin_interp.idx_function_real_grid)
        
        # Linearly Interpolate between splines
        return LinInterp.interpolate(x_imag, compl_spline_lin_interp.imaginaries[idx], compl_spline_lin_interp.imaginaries[idx+1], f_real_prev, f_real_next) +
               1im * LinInterp.interpolate(x_imag, compl_spline_lin_interp.imaginaries[idx], compl_spline_lin_interp.imaginaries[idx+1], f_imag_prev, f_imag_next)
    end
end