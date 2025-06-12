module CSplineInterpolation
    using LinearAlgebra
    
    struct CSplineInterpolator
        x::AbstractArray{<:Number}

        y::AbstractArray{<:Number}
        dd_y::AbstractArray{<:Number}

        CSplineInterpolator(x::AbstractArray{<:Number}, f::AbstractArray{<:Number}) = begin
            dd_f = solve_interpolation(x, f)
            return new(x, f, dd_f)
        end
    end

    function solve_interpolation(x::AbstractArray{<:Number}, y::AbstractArray{<:Number})::AbstractArray{Number}
        N = length(x)
        
        r = zeros(Float64, N)
        a = zeros(Float64, N)
        b = zeros(Float64, N)
        c = zeros(Float64, N)

        for i = 2:(N-1)
            r[i] = (y[i+1]-y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1])
            a[i] = (x[i] - x[i-1]) / 6.0
            b[i] = (x[i+1] - x[i-1]) / 3.0
            c[i] = (x[i+1] - x[i]) / 6.0
        end

        r[1] = 0
        b[1] = 1
        c[1] = 0

        r[N] = 0
        a[N] = 0
        b[N] = 1


        # Solve Matrix Equation by forward and backward substitution
        beta = zeros(Float64, N)
        rho = zeros(Float64, N)
        
        beta[1] = b[1]
        rho[1] = r[1]
        for i = 2:N
            beta[i] = b[i] - a[i] / beta[i-1] * c[i-1]
            rho[i] = r[i] - a[i] / beta[i-1] * rho[i-1]
        end

        z = zeros(Float64, N)
        z[N] = rho[N] / beta[N]

        for j = 1:N-1
            z[N-j] = (rho[N-j] - c[N-j] * z[N-j+1]) / beta[N-j]
        end

        return z
    end

    # Spacing function yields how the grid is interpolated, it should have signature
    #   idx_function(x, num_points, lower_bound, upper_bound) and return an index i s.t. (x_i < x < x_i+1)
    function interpolate(x::Number, cspline_interp::CSplineInterpolator, idx_function::Function)

        # TODO move if to "parent" method (same for Linear Interpolation)
        if (cspline_interp.x[1] < x && x < cspline_interp.x[end])
            i::Int = idx_function(x, cspline_interp.x)

            if (i < 0)
                print(idx_function)
            end
        
            xi = cspline_interp.x[i]
            xip1 = cspline_interp.x[i+1]

            A_cur = A(x, xi, xip1)
            B_cur = B(x, xi, xip1)

            return A_cur * cspline_interp.y[i] + B_cur * cspline_interp.y[i+1] + CD(A_cur, xi, xip1) * cspline_interp.dd_y[i] + CD(B_cur, xi, xip1) * cspline_interp.dd_y[i+1]
            
        elseif (x == cspline_interp.x[1])
            return cspline_interp.y[1]

        elseif (x == cspline_interp.x[end])
            return cspline_interp.y[end]
        end

        throw(ErrorException("Out of bounds for interpolation: $(x) not in [$(maximum(cspline_interp.x)), $(minimum(cspline_interp.x))]"))
        
    end

    function CD(AB::Float64, xi::Float64, xip1::Float64)
        return 1.0/6.0 * (AB^3 - AB) * (xip1 - xi)^2
    end

    function A(x, xi, xip1)
        return (xip1 - x) / (xip1 - xi)
    end

    function B(x, xi, xip1)
        return (x - xi) / (xip1 - xi)
    end
end