module GaussLegendre

    import ..Integrators
    import ....Numerics.LegendrePolynomials as Polynomials
    import ....Numerics.NewtonRaphson as NewtonRaphson

    ##################### Gauss Legendre #####################
    struct GaussLegendreData
        x::AbstractArray{Float64}
        w::AbstractArray{Float64}

        GaussLegendreData(N::Int) = begin
            narray = Float64.(1:N)
            x = zeros(N)

            # generate legpol
            leg_dleg_N(x) = Polynomials.legendre_poly_and_deriv(x, N)

            # find zeros
            x_guesses = cos.((narray .- 0.25) .* (pi / (N + 0.5)))
            
            len = length(x_guesses)
            for i = 1:len
                x[len-i+1] = NewtonRaphson.find_roots(leg_dleg_N, x_guesses[i], 1E-12)
            end

            # generate weights
            w = 2.0 / (N + 1.0)^2 .* ((1.0 .- x.^2) ./ (Polynomials.legendre_poly.(x, N+1)).^2)

            return new(x, w)
        end
    end    

    struct GaussLegendreIntegrator <: Integrators.Integrator
        data::GaussLegendreData

        GaussLegendreIntegrator(N::Int) = begin
            return new(GaussLegendreData(N))
        end
    end

    function integrate(integrator::GaussLegendreIntegrator, f::Function)
        sum::ComplexF64 = 0
        for i = 1:length(integrator.data.x)
            sum += integrator.data.w[i] * f(integrator.data.x[i])
        end

        return sum
    end

    function integrate(integrator::GaussLegendreIntegrator, f::Function, a::Number, b::Number)
        sum::ComplexF64 = 0
        for i = 1:length(integrator.data.x)
            sum += integrator.data.w[i] * f(((b - a) * integrator.data.x[i] + a + b) / 2.0)
        end

        return (b - a) / 2.0 * sum
    end

    function integrate_logspaced(integrator::GaussLegendreIntegrator, f::Function, a::Number, b::Number)
        sum::ComplexF64 = 0.0
        
        A = - log(a * b) / log(b / a)
        B = 2.0 / log(b / a)

        # TODO
        for i = 1:length(integrator.data.x)
            z = exp((integrator.data.x[i] - A)/B)
            sum += f(z) * z * integrator.data.w[i]
        end

        return sum / B
    end

    function get_logspaced_grid_matching(integrator::GaussLegendreIntegrator, a::Number, b::Number)  # TODO move to base module
        A = - log(a * b) / log(b / a)
        B = 2.0 / log(b / a)
        
        grid = zeros(typeof(integrator.data.x[1]), size(integrator.data.x))
        for i in eachindex(integrator.data.x)
            grid[i] = exp((integrator.data.x[i] - A) / B)
        end

        return grid
    end
end