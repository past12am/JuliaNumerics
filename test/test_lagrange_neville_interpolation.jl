# =============================================================================
# Validation of LagrangeNevilleInterpolation (Neville's algorithm, NR3 §3.2)
# -----------------------------------------------------------------------------
# The interpolant through M nodes is the unique polynomial of degree M - 1, so
# it must reproduce such polynomials exactly, hit the data at the nodes, agree
# with the explicit Lagrange form, and not depend on the node order. A 4-point
# stencil must converge as h⁴.
#
# Run from the JuliaNumerics directory with
#   julia test/test_lagrange_neville_interpolation.jl [<path to Numerics/src>]
# =============================================================================

using Test

const NUMERICS_SRC = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "src")
include(joinpath(NUMERICS_SRC, "Numerics.jl"))

const LN = Numerics.Interpolators.LagrangeNevilleInterpolation


# Independent reference: explicit Lagrange form
function lagrange_reference(x_vals, y_vals, x)
    s = zero(promote_type(eltype(y_vals), typeof(x)))
    for j in eachindex(x_vals)
        l = one(x)
        for k in eachindex(x_vals)
            k == j && continue
            l *= (x - x_vals[k]) / (x_vals[j] - x_vals[k])
        end
        s += l * y_vals[j]
    end
    return s
end

# Non-uniform, deterministic nodes on roughly [0, 2]
nodes(M) = [0.3 * k + 0.05 * sin(3.0 * k) for k = 0:M-1]

f_test(x) = exp(x) * sin(3.0 * x)

# Evaluation points inside the stencil and slightly beyond it (extrapolation)
eval_points(x_vals) = range(x_vals[1] - 0.2, x_vals[end] + 0.2, length=41)


@testset "LagrangeNevilleInterpolation" begin

    @testset "single node" begin
        @test LN.interpolate([0.7], [2.5], 1.3) == 2.5
    end

    @testset "polynomials of degree M - 1 are exact" begin
        coeffs = [0.4, -1.3, 0.8, 2.1, -0.6, 0.25, -0.1]
        for M = 1:7
            p(x) = sum(coeffs[j+1] * x^j for j = 0:M-1)
            x_vals = nodes(M)
            y_vals = p.(x_vals)
            for x in eval_points(x_vals)
                @test LN.interpolate(x_vals, y_vals, x) ≈ p(x) rtol = 1e-12 atol = 1e-12
            end
        end
    end

    @testset "nodes are reproduced" begin
        for M = 2:7
            x_vals = nodes(M)
            y_vals = f_test.(x_vals)
            for k = 1:M
                @test LN.interpolate(x_vals, y_vals, x_vals[k]) ≈ y_vals[k] rtol = 1e-13 atol = 1e-14
            end
        end
    end

    @testset "agrees with the Lagrange form" begin
        for M = 2:7
            x_vals = nodes(M)
            y_vals = f_test.(x_vals)
            for x in eval_points(x_vals)
                @test LN.interpolate(x_vals, y_vals, x) ≈ lagrange_reference(x_vals, y_vals, x) rtol = 1e-12 atol = 1e-13
            end
        end
    end

    @testset "node order is irrelevant" begin
        x_vals = nodes(5)
        y_vals = f_test.(x_vals)
        perm = [3, 5, 1, 4, 2]
        for x in eval_points(x_vals)
            @test LN.interpolate(x_vals[perm], y_vals[perm], x) ≈ LN.interpolate(x_vals, y_vals, x) rtol = 1e-12 atol = 1e-13
        end
    end

    @testset "complex values, real nodes" begin
        x_vals = nodes(4)
        y_vals = ComplexF64.(f_test.(x_vals), cos.(x_vals))
        for x in eval_points(x_vals)
            y_re = LN.interpolate(x_vals, real.(y_vals), x)
            y_im = LN.interpolate(x_vals, imag.(y_vals), x)
            y = LN.interpolate(x_vals, y_vals, x)
            @test y isa ComplexF64
            @test y ≈ complex(y_re, y_im) rtol = 1e-13 atol = 1e-14
        end
    end

    @testset "4-point stencil on a view of a log grid" begin
        # Intended use: 4 neighbouring nodes in ln p⃗², passed as views
        lnp2_grid = collect(range(log(1e-3), log(1e3), length=64))
        g_grid = f_test.(lnp2_grid ./ 4)
        for j = 1:length(lnp2_grid)-3
            x_stencil = @view lnp2_grid[j:j+3]
            y_stencil = @view g_grid[j:j+3]
            x = 0.5 * (x_stencil[2] + x_stencil[3])
            @test LN.interpolate(x_stencil, y_stencil, x) ≈ lagrange_reference(x_stencil, y_stencil, x) rtol = 1e-12 atol = 1e-13
        end
    end

    @testset "4-point stencil converges as h⁴" begin
        # Midpoint of the centre interval: error = f⁗(ξ)/4! · (9/16) h⁴, f⁗ ≠ 0 near x0
        x0 = 1.0
        err(h) = begin
            x_vals = [x0 + h * k for k = -1:2]
            abs(LN.interpolate(x_vals, f_test.(x_vals), x0 + h / 2) - f_test(x0 + h / 2))
        end
        for h in (0.1, 0.05, 0.025)
            @test 14.0 < err(h) / err(h / 2) < 18.0
        end
    end

    @testset "interpolate! with preallocated buffers" begin
        x_vals = nodes(4)
        y_real = f_test.(x_vals)
        y_complex = ComplexF64.(y_real, cos.(x_vals))

        C = zeros(ComplexF64, 4)
        D = zeros(ComplexF64, 4)

        # Reused buffers: real and complex data alternately, no stale state
        for x in eval_points(x_vals)
            @test LN.interpolate!(C, D, x_vals, y_complex, x) ≈ LN.interpolate(x_vals, y_complex, x) rtol = 1e-14 atol = 1e-15
            @test LN.interpolate!(C, D, x_vals, y_real, x) ≈ LN.interpolate(x_vals, y_real, x) rtol = 1e-14 atol = 1e-15
        end

        # Allocation-free on views of a larger grid (warm-up call compiles first)
        lnp2_grid = collect(range(log(1e-3), log(1e3), length=64))
        g_grid = ComplexF64.(f_test.(lnp2_grid ./ 4), cos.(lnp2_grid))
        x_stencil = @view lnp2_grid[20:23]
        y_stencil = @view g_grid[20:23]
        x = 0.5 * (x_stencil[2] + x_stencil[3])

        count_allocations(C, D, xs, ys, x) = @allocated LN.interpolate!(C, D, xs, ys, x)
        count_allocations(C, D, x_stencil, y_stencil, x)
        @test count_allocations(C, D, x_stencil, y_stencil, x) == 0
    end

    @testset "identical nodes throw" begin
        @test_throws ArgumentError LN.interpolate([0.1, 0.5, 0.5, 0.9], [1.0, 2.0, 3.0, 4.0], 0.3)
    end
end
