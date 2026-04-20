module Numerics
    include("PolynomialsAndFunctions/LegendrePolynomials.jl")
    include("PolynomialsAndFunctions/ChebyshevPolynomials.jl")
    include("PolynomialsAndFunctions/ModifiedBesselFirstKind.jl")

    include("Roots/NewtonRaphson.jl")
    include("Integration/Integrators.jl")

    include("Interpolation/Interpolators.jl")

    include("Extrapolation/Schlessinger.jl")
end