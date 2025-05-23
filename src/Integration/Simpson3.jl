module Simpson3
    import ..Integrators
    
    ##################### Simpson 3 Integration #####################
    struct SimpsonIntegrator <: Integrators.Integrator 
        N::Int

        SimpsonIntegrator(N::Int) = begin
            if N % 2 == 0
                throw(ErrorException("Need an odd number of grid-points for Simpson Integrator"))
            end

            return new(N)
        end
    end

    function integrate(integrator::SimpsonIntegrator, f::Function, a::Float64, b::Float64)
        res::Float64 = f(a) + f(b)

        h::Float64 = (b - a) / (integrator.N - 1)
        h_twice::Float64 = 2.0 * h

        c1::Float64 = 0.0
        x_cur = a + h
        for i = 2:2:(integrator.N-1)
            c1 += f(x_cur)
            x_cur += h_twice
        end

        c2::Float64 = 0.0
        x_cur = a + h_twice
        for i = 3:2:(integrator.N-3)
            c2 += f(x_cur)
            x_cur += h_twice
        end

        res += 4.0 * c1 + 2.0 * c2
        res *= h / 3.0
        
        return res
    end
end