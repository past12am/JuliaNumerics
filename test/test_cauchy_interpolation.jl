# =============================================================================
# Validation of CauchyInterpolatorParabolaSchwarzReflectable
# -----------------------------------------------------------------------------
# Contour: upper parabola z(t) = t - M²/4 + i M √t, t ∈ [0, Λ²] (Gauss-Legendre
# linear in q = √t on [0, √switch2], log-spaced in t on [switch2, Λ²]), its
# Schwarz reflection (lower branch), and the vertical connection at
# Re z = Λ² - M²/4.
#
# Test functions are analytic inside and on the contour and satisfy
# f(z̄) = conj f(z), so Cauchy's formula must reproduce them inside. The
# barycentric form is only valid inside the contour; outside it must throw.
#
# Run from the JuliaNumerics directory with
#   julia test/test_cauchy_interpolation.jl [<path to Numerics/src>]
# =============================================================================

using Test

const NUMERICS_SRC = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "src")
include(joinpath(NUMERICS_SRC, "Numerics.jl"))

const CI = Numerics.Interpolators.CauchyIntegralInterpolation


# Point inside the parabola for |y| < 1 (y = ±1 is the contour)
z_of_qy(q, y, M) = ComplexF64(q^2 - M^2 / 4, M * q * y)

# Strictly inside the closed contour
inside(z, M, Λ2) = real(z) < Λ2 - M^2 / 4 && imag(z)^2 < M^2 * (real(z) + M^2 / 4)

function cauchy(f, z, ip)
    f_grid = ComplexF64.(f.(ip.z_at_gridpoint))
    f_inf = ComplexF64.(f.(ip.z_at_infconnection))
    return CI.interpolate(ComplexF64(z), f_grid, f_inf, ip)
end

# Pole pair at z0, z̄0 must lie outside the parabola for every M used below
const z0 = -0.6 + 0.8im
const test_functions = [
    ("constant",     z -> 1.0 + 0.0im),
    ("polynomial",   z -> z^2 + 0.3z - 1.0),
    ("real pole",    z -> 1.0 / (z + 1.0)),
    ("pole pair",    z -> 1.0 + 0.5 / ((z - z0) * (z - conj(z0)))),
    ("logarithm",    z -> log(z + 2.0)),
]

const Λ2 = 1e3
const switch2 = 1.0
const N_lin = 64
const N_log = 192
const N_parabola = N_lin + N_log
const N_inf = 64

# Plain (non-barycentric) Cauchy is off by O(1e-1 … 1) at |y| = 0.99
const tol_inside = 1e-7

# Up to the contour (y = ±1) and into the UV; the apex region is tested separately
const qs = (0.3, 1.0, 3.0, 10.0, 25.0)
const ys = (-0.99, -0.9, -0.5, 0.0, 0.5, 0.9, 0.99)


@testset "CauchyInterpolatorParabolaSchwarzReflectable" begin

    for M in (0.5, 1.0, 1.6)
        ip = CI.CauchyInterpolatorParabolaSchwarzReflectable(N_lin, N_log, N_inf, switch2, Λ2, M)

        @testset "geometry (M = $M)" begin
            # Nodes lie on the upper parabola z(t) = t - M²/4 + i M √t, t ∈ (0, Λ²), ascending
            @test all(0.0 .< ip.t_at_gridpoint .< Λ2)
            @test issorted(ip.t_at_gridpoint)
            @test count(ip.t_at_gridpoint .< switch2) == N_lin
            for (t, z) in zip(ip.t_at_gridpoint, ip.z_at_gridpoint)
                @test real(z) + M^2 / 4 ≈ t atol = 1e-14 * M^2
                @test imag(z) ≈ M * sqrt(t) rtol = 1e-14
            end
            # Contour weights integrate dz from the apex: Σ w_i = z(Λ²) - z(0), on each piece too
            z_end = CI.z_of_t(Λ2, M)
            @test sum(ip.w_at_gridpoint) ≈ z_end - CI.z_of_t(0.0, M) rtol = 1e-12
            @test sum(ip.w_at_gridpoint[1:N_lin]) ≈ CI.z_of_t(switch2, M) - CI.z_of_t(0.0, M) rtol = 1e-12
            # Connection closes the contour at the parabola's end point z(Λ²), upwards
            @test all(real.(ip.z_at_infconnection) .≈ real(z_end))
            @test all(abs.(imag.(ip.z_at_infconnection)) .< imag(z_end))
            @test minimum(imag.(ip.z_at_infconnection)) < 0 < maximum(imag.(ip.z_at_infconnection))
            @test sum(ip.w_at_infconnection) ≈ 2im * imag(z_end) rtol = 1e-12
        end

        @testset "reproduces analytic functions inside (M = $M)" begin
            for (name, f) in test_functions, q in qs, y in ys
                z = z_of_qy(q, y, M)
                @test inside(z, M, Λ2)
                @test isapprox(cauchy(f, z, ip), f(z); rtol = tol_inside, atol = tol_inside)
            end
        end

        @testset "Schwarz reflection (M = $M)" begin
            for (name, f) in test_functions, q in qs, y in (0.5, 0.99)
                z = z_of_qy(q, y, M)
                @test cauchy(f, conj(z), ip) ≈ conj(cauchy(f, z, ip)) rtol = 1e-9
            end
            # Real axis inside the parabola, including (-M²/4, 0): result is real
            for x in (-0.5 * M^2 / 4, 0.0, 0.5, 5.0, 50.0, 500.0), (name, f) in test_functions
                val = cauchy(f, x, ip)
                # exact up to the symmetry of the GL nodes on the connection
                @test abs(imag(val)) <= 1e-8 * max(1.0, abs(val))
                @test isapprox(val, f(x); rtol = tol_inside, atol = tol_inside)
            end
        end

        @testset "value at the nodes (M = $M)" begin
            # Barycentric form is 0/0 exactly on a node: must return the stored value
            for (name, f) in test_functions, i in (1, N_parabola ÷ 2, N_parabola)
                f_grid = ComplexF64.(f.(ip.z_at_gridpoint))
                f_inf = ComplexF64.(f.(ip.z_at_infconnection))
                @test CI.interpolate(ip.z_at_gridpoint[i], f_grid, f_inf, ip) == f_grid[i]
                @test CI.interpolate(conj(ip.z_at_gridpoint[i]), f_grid, f_inf, ip) == conj(f_grid[i])
            end
            # Points on the connection are inside the check's tolerance: same on its nodes
            for (name, f) in test_functions, i in (1, N_inf ÷ 2, N_inf)
                f_grid = ComplexF64.(f.(ip.z_at_gridpoint))
                f_inf = ComplexF64.(f.(ip.z_at_infconnection))
                @test CI.interpolate(ip.z_at_infconnection[i], f_grid, f_inf, ip) == f_inf[i]
            end
        end

        @testset "several functions in one call (M = $M)" begin
            # Tuple version: shared r_i, each entry bit-identical to the single-function call
            f_grid = Tuple(ComplexF64.(f.(ip.z_at_gridpoint)) for (name, f) in test_functions)
            f_inf = Tuple(ComplexF64.(f.(ip.z_at_infconnection)) for (name, f) in test_functions)
            for q in qs, y in ys
                z = z_of_qy(q, y, M)
                vals = CI.interpolate(z, f_grid, f_inf, ip)
                @test vals isa NTuple{length(test_functions), ComplexF64}
                @test all(vals[k] === CI.interpolate(z, f_grid[k], f_inf[k], ip) for k in eachindex(test_functions))
            end
            # Nodes (parabola, its reflection, connection) return the stored values
            for i in (1, N_parabola)
                @test CI.interpolate(ip.z_at_gridpoint[i], f_grid, f_inf, ip) == getindex.(f_grid, i)
                @test CI.interpolate(conj(ip.z_at_gridpoint[i]), f_grid, f_inf, ip) == conj.(getindex.(f_grid, i))
            end
            for i in (1, N_inf ÷ 2, N_inf)
                @test CI.interpolate(ip.z_at_infconnection[i], f_grid, f_inf, ip) == getindex.(f_inf, i)
            end
        end

        @testset "outside the contour throws (M = $M)" begin
            outside = (ComplexF64(-M^2, 0.0),                 # left of the apex
                       ComplexF64(2 * Λ2, 0.0),               # beyond the connection
                       ComplexF64(1.0, 1.01 * M * sqrt(1.0 + M^2 / 4)),
                       ComplexF64(1.0, -1.01 * M * sqrt(1.0 + M^2 / 4)))
            for z in outside
                @test !inside(z, M, Λ2)
                @test_throws DomainError cauchy(z -> 1.0 + 0.0im, z, ip)
            end
        end

        @testset "apex region (M = $M)" begin
            # Points very close to the apex (dist ≳ 1e-6), where the contour is parametrized by q
            for (name, f) in test_functions
                for z in (ComplexF64(-0.9 * M^2 / 4, 0.0), z_of_qy(0.01, 0.0, M), z_of_qy(0.01, 0.9, M),
                          z_of_qy(0.05, 0.0, M), z_of_qy(0.05, 0.99, M))
                    @test isapprox(cauchy(f, z, ip), f(z); rtol = 1e-8, atol = 1e-8)
                end
            end
        end
    end

    @testset "convergence in the number of parabola points" begin
        M = 1.0
        f = z -> 1.0 / (z + 1.0)
        z = z_of_qy(1.0, 0.9, M)
        errs = [abs(cauchy(f, z, CI.CauchyInterpolatorParabolaSchwarzReflectable(N ÷ 4, N - N ÷ 4, N_inf, switch2, Λ2, M)) - f(z))
                for N in (16, 32, 64)]
        @test errs[2] < errs[1]
        @test errs[3] < errs[2]
    end
end
