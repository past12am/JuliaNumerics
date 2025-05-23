
module Integrators
    abstract type Integrator end
    integrate(::Integrator, f::Function) = error("integrate not implemented for $(type(Integrator))")
    integrate(::Integrator, f::Function, a::Float64, b::Float64) = error("integrate not implemented for $(type(Integrator))")
    
    include("GaussChebyshev.jl")
    include("GaussLegendre.jl")
    include("Simpson3.jl")
end