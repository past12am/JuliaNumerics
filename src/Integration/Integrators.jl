
module Integrators
    abstract type Integrator end
    integrate(integrator::Integrator, f::Function) = error("integrate not implemented for $(typeof(integrator))")
    integrate(integrator::Integrator, f::Function, a::Float64, b::Float64) = error("integrate not implemented for $(typeof(integrator))")
    
    include("GaussChebyshev.jl")
    include("GaussLegendre.jl")
    include("Simpson3.jl")
end