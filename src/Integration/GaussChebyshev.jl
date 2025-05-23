module GaussChebyshev
    import ..Integrators
    
    ##################### Gauss Chebyshev #####################
    struct GaussChebyData
        x::AbstractArray{Float64}
        w::AbstractArray{Float64}

        GaussChebyData(N::Int) = begin
            narray = Float64.(1:N)

            x = cos.(narray .* (pi / (N + 1.0)))
            w = pi / (N + 1) .* sin.(narray .* (pi / (N + 1.0))).^2

            return new(x, w)
        end
    end    

    struct GaussChebyIntegrator <: Integrators.Integrator
        data::GaussChebyData

        GaussChebyIntegrator(N::Int) = begin
            return new(GaussChebyData(N))
        end
    end

    function integrate(integrator::GaussChebyIntegrator, f::Function)
        res::Float64 = 0.0
        for i = 1:length(integrator.data.x)
            res += integrator.data.w[i] * f(integrator.data.x[i])
        end

        return res
    end
end