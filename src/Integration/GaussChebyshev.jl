module GaussChebyshev
    import ..Integrators
    
    ##################### Gauss Chebyshev #####################
    struct GaussChebyData
        x::AbstractArray{Float64}
        w::AbstractArray{Float64}

        GaussChebyData(N::Int) = begin
            narray = Float64.(1:N)

            a = pi / (N + 1.0)

            x = cos.(narray .* a)
            w = a .* (sin.(narray .* a)).^2

            return new(x, w)
        end
    end    

    struct GaussChebyIntegrator <: Integrators.Integrator
        data::GaussChebyData

        GaussChebyIntegrator(N::Int) = begin
            return new(GaussChebyData(N))
        end
    end

    # TODO check GC integrator
    function integrate(integrator::GaussChebyIntegrator, f::Function)
        res::Float64 = 0.0
        for i = 1:length(integrator.data.x)
            res += integrator.data.w[i] * f(integrator.data.x[i])
        end

        return res
    end

    function integrate_complex(integrator::GaussChebyIntegrator, f::Function)
        res::ComplexF64 = 0.0
        for i = 1:length(integrator.data.x)
            res += integrator.data.w[i] * f(integrator.data.x[i])
        end

        return res
    end
end