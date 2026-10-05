module CauchyIntegralInterpolation
    import ....Numerics.Integrators.GaussLegendre as GL
    

    function z_of_t(t::Float64, mass_limit::Float64)::ComplexF64
        return ComplexF64(t - mass_limit*mass_limit / 4.0, sqrt(t) * mass_limit)
    end

    function dz_dt_of_t(t::Float64, mass_limit::Float64)::ComplexF64
        ComplexF64(1.0, mass_limit / (2.0 * sqrt(t)))
    end


    function z_of_t_infconnection(t::Float64, re_part::Float64)::ComplexF64
        return ComplexF64(re_part, t)
    end


    struct CauchyInterpolatorParabolaSchwarzReflectable
    
        mass_limit2::Float64
        re_part_infconnection::Float64

        Λ2::Float64
        
        gl_parabola_integrator_ir::GL.GaussLegendreIntegrator
        gl_parabola_integrator_mid_to_uv::GL.GaussLegendreIntegrator
        
        t_at_gridpoint::Vector{Float64}
        z_at_gridpoint::Vector{ComplexF64}
        w_at_gridpoint::Vector{ComplexF64}

        
        gl_infconnection_integrator::GL.GaussLegendreIntegrator

        t_at_infconnection::Vector{Float64}
        z_at_infconnection::Vector{ComplexF64}
        w_at_infconnection::Vector{ComplexF64}

        CauchyInterpolatorParabolaSchwarzReflectable(num_parabola_points_linspaced::Int, num_parabola_points_logspaced::Int, num_infconnection_points::Int, lin_to_log_switch2::Float64, Λ2::Float64, mass_limit::Float64) = begin
            mass_limit2::Float64 = mass_limit*mass_limit
            
            re_part_infconnection::Float64 = Λ2 - mass_limit2 / 4.0

            # Linear Grid in the IR     (linear in sqrt(t), not t)
            gl_parabola_integrator_ir::GL.GaussLegendreIntegrator = GL.GaussLegendreIntegrator(num_parabola_points_linspaced)

            sqrt_t_at_gridpoint_ir::Vector{Float64} = GL.get_linear_grid_matching(gl_parabola_integrator_ir, 0.0, sqrt(lin_to_log_switch2))
            sqrt_t_jacobian_at_gridpoint_ir::Vector{Float64} = GL.get_linear_jacobian_matching(gl_parabola_integrator_ir, 0.0, sqrt(lin_to_log_switch2))
            
            t_at_gridpoint_ir::Vector{Float64} = sqrt_t_at_gridpoint_ir.^2
            jacobian_at_gridpoint_ir::Vector{Float64} = 2.0 .* sqrt_t_at_gridpoint_ir .* sqrt_t_jacobian_at_gridpoint_ir    # dt = 2 q dq   with q = sqrt(t)
            
            z_at_gridpoint_ir::Vector{ComplexF64} = z_of_t.(t_at_gridpoint_ir, mass_limit)
            dz_dt_at_gridpoint_ir::Vector{ComplexF64} = dz_dt_of_t.(t_at_gridpoint_ir, mass_limit)
            
            w_at_gridpoint_ir::Vector{ComplexF64} = gl_parabola_integrator_ir.data.w .* jacobian_at_gridpoint_ir .* dz_dt_at_gridpoint_ir


            # Logspaced grid from MID to UV
            gl_parabola_integrator_mid_to_uv::GL.GaussLegendreIntegrator = GL.GaussLegendreIntegrator(num_parabola_points_logspaced)

            t_at_gridpoint_mid_to_uv::Vector{Float64} = GL.get_logspaced_grid_matching(gl_parabola_integrator_mid_to_uv, lin_to_log_switch2, Λ2)
            jacobian_at_gridpoint_mid_to_uv::Vector{Float64} = GL.get_logspaced_jacobian_matching(gl_parabola_integrator_mid_to_uv, lin_to_log_switch2, Λ2)

            z_at_gridpoint_mid_to_uv::Vector{ComplexF64} = z_of_t.(t_at_gridpoint_mid_to_uv, mass_limit)
            dz_dt_at_gridpoint_mid_to_uv::Vector{ComplexF64} = dz_dt_of_t.(t_at_gridpoint_mid_to_uv, mass_limit)

            w_at_gridpoint_mid_to_uv::Vector{ComplexF64} = gl_parabola_integrator_mid_to_uv.data.w .* jacobian_at_gridpoint_mid_to_uv .* dz_dt_at_gridpoint_mid_to_uv
            

            t_at_gridpoint = vcat(t_at_gridpoint_ir, t_at_gridpoint_mid_to_uv)
            z_at_gridpoint = vcat(z_at_gridpoint_ir, z_at_gridpoint_mid_to_uv)
            w_at_gridpoint = vcat(w_at_gridpoint_ir, w_at_gridpoint_mid_to_uv)


            # Connection at Infinity = cutoff2
            gl_infconnection_integrator::GL.GaussLegendreIntegrator = GL.GaussLegendreIntegrator(num_infconnection_points)

            t_at_infconnection::Vector{Float64} = GL.get_linear_grid_matching(gl_infconnection_integrator, -sqrt(Λ2) * mass_limit, sqrt(Λ2) * mass_limit)
            jacobian_at_infconnection::Vector{Float64} = GL.get_linear_jacobian_matching(gl_infconnection_integrator, -sqrt(Λ2) * mass_limit, sqrt(Λ2) * mass_limit)

            z_at_infconnection::Vector{ComplexF64} = z_of_t_infconnection.(t_at_infconnection, re_part_infconnection)
            dz_dt_at_infconnection = ComplexF64(0.0, 1.0)

            w_at_infconnection::Vector{ComplexF64} = gl_infconnection_integrator.data.w .* jacobian_at_infconnection .* dz_dt_at_infconnection

            new(mass_limit2, re_part_infconnection, Λ2, gl_parabola_integrator_ir, gl_parabola_integrator_mid_to_uv, t_at_gridpoint, z_at_gridpoint, w_at_gridpoint, gl_infconnection_integrator, t_at_infconnection, z_at_infconnection, w_at_infconnection)
        end
    end


    # w / d without the overflow/underflow scaling of Base's complex division (|d| is bounded on and inside the contour)
    @inline function div_unscaled(w::ComplexF64, d::ComplexF64)::ComplexF64
        return w * conj(d) * (1.0 / abs2(d))
    end


    function interpolate(a::ComplexF64, f_at_gridpoint::NTuple{K, Vector{ComplexF64}}, f_at_infconnection::NTuple{K, Vector{ComplexF64}}, interp::CauchyInterpolatorParabolaSchwarzReflectable)::NTuple{K, ComplexF64} where {K}
        
        # check that a lies inside the contour
        re_shifted::Float64 = a.re + interp.mass_limit2 / 4.0
        tol::Float64 = 1e-12 * (abs(a.re) + interp.mass_limit2)
        if (a.re > interp.Λ2 - interp.mass_limit2 / 4.0 + tol || a.im * a.im > interp.mass_limit2 * (re_shifted + tol))
            throw(DomainError(a, "outside the interpolation parabola"))
        end

        a_conj::ComplexF64 = conj(a)

        # Integrate logspaced along contour
        gamma_1_contour_res::NTuple{K, ComplexF64} = ntuple(_ -> zero(ComplexF64), Val(K))
        gamma_2_contour_res::NTuple{K, ComplexF64} = ntuple(_ -> zero(ComplexF64), Val(K))

        gamma_1_simple_pole_integ::ComplexF64 = zero(ComplexF64)
        gamma_2_simple_pole_integ::ComplexF64 = zero(ComplexF64)

        for i = 1:length(interp.t_at_gridpoint)
            z_i = interp.z_at_gridpoint[i]

            # If we sit on the node, return the value directly
            if z_i == a
                return ntuple(k -> f_at_gridpoint[k][i], Val(K))
            elseif z_i == a_conj
                return ntuple(k -> conj(f_at_gridpoint[k][i]), Val(K))
            end

            w_i = interp.w_at_gridpoint[i]
            r_i = div_unscaled(w_i, (z_i - a))
            r_i_conj = div_unscaled(w_i, (z_i - a_conj))

            f_i = ntuple(k -> f_at_gridpoint[k][i], Val(K))

            gamma_1_contour_res = gamma_1_contour_res .+ (r_i .* f_i)
            gamma_2_contour_res = gamma_2_contour_res .+ (r_i_conj .* f_i)

            gamma_1_simple_pole_integ += r_i
            gamma_2_simple_pole_integ += r_i_conj
        end

        #   Fix contour direction (counter-clockwise) and conjugation
        gamma_1_contour_res = -1.0 .* gamma_1_contour_res
        gamma_2_contour_res = conj.(gamma_2_contour_res)

        gamma_1_simple_pole_integ *= -1.0
        gamma_2_simple_pole_integ = conj(gamma_2_simple_pole_integ)

        


        # Integrate Linear at infconnection
        gamma_3_contour_res::NTuple{K, ComplexF64} = ntuple(_ -> zero(ComplexF64), Val(K))

        gamma_3_simple_pole_integ::ComplexF64 = zero(ComplexF64)

        for i = 1:length(interp.t_at_infconnection)
            z_i = interp.z_at_infconnection[i]

            # If we sit on the node, return the value directly
            if z_i == a
                return ntuple(k -> f_at_infconnection[k][i], Val(K))
            end

            r_i = div_unscaled(interp.w_at_infconnection[i], (z_i - a))
            f_i = ntuple(k -> f_at_infconnection[k][i], Val(K))

            gamma_3_contour_res = gamma_3_contour_res .+ (r_i .* f_i)

            gamma_3_simple_pole_integ += r_i
        end

        contour_res = gamma_1_contour_res .+ gamma_2_contour_res .+ gamma_3_contour_res
        simple_pole_res = gamma_1_simple_pole_integ + gamma_2_simple_pole_integ + gamma_3_simple_pole_integ

        # Cauchy Integral Formula to interpolate (note counter-clockwise orientation of contour)
        return contour_res ./ simple_pole_res
    end

    function interpolate(a::ComplexF64, f_at_gridpoint::Vector{ComplexF64}, f_at_infconnection::Vector{ComplexF64}, interp::CauchyInterpolatorParabolaSchwarzReflectable)::ComplexF64
        return interpolate(a, (f_at_gridpoint,), (f_at_infconnection,), interp)[1]
    end
end