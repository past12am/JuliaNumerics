module ModifiedBesselfunctionsFirstKind

    # Assuming n_request_list is sorted
    function bessel_I_n(x::Number, n_request_list::AbstractArray{Int})
        res = zeros(size(n_request_list))
        handled_n = falses(size(n_request_list))

        greater2idxs = (n_request_list .>= 2)
        handled_n .= .!greater2idxs

        if any(greater2idxs)
            
            # Asymtotics
            for (idx, n) in enumerate(n_request_list)
                if (n < 2)
                    continue
                end

                if x < sqrt(n) / 2          # TODO better condition
                    handled_n[idx] = true
                    res[idx] = bessel_I_n_asympt_xsmall(x, n)
                elseif x > 2 * n            # TODO better condition
                    handled_n[idx] = true
                    res[idx] = bessel_I_n_asympt_xlarge(x)
                end
            end

            # Handle remaining
            if !all(handled_n)
                acc::Float64 = 200.0
                n_max = maximum(n_request_list)
                n_start::Int = floor(Int, 2 * ((n_max + 1) + sqrt(acc * (n_max + 1))))

                bessel_p::Float64 = 0.0
                bessel::Float64 = 1.0
                bessel_m::Float64 = 0.0

                ctr::Int = 0

                two_over_x = 2.0 / abs(x)
                for n = n_start:-1:1
                    bessel_m = bessel_p + n * two_over_x * bessel

                    bessel_p = bessel
                    bessel = bessel_m

                    if (bessel > 1E10)
                        bessel *= 1E-10
                        bessel_p *= 1E-10
                        bessel_m *= 1E-10

                        res .*= 1E-10
                    end

                    if n >= 2 && ctr < length(res) && n == (n_request_list[end - ctr] + 1) && ! handled_n[end - ctr]
                        res[end - ctr] = bessel_m  # we already shifted
                        ctr += 1
                    end
                end

                scaling = i0(x) / bessel_m
                res .*= scaling
            end
        end

        # For the lowest 2, use our special functions
        zero_idxs = findall((n) -> n == 0, n_request_list)
        one_idxs = findall((n) -> n == 1, n_request_list)

        res[zero_idxs] .= i0(x)
        res[one_idxs] .= i1(x)

        return res
    end


    function poly(cof::Vector{Float64}, n::Int, x::Number)
        ans = cof[n+1]
        for i in n:-1:1
            ans = ans * x + cof[i]
        end
        return ans
    end

    i0p = [9.999999999999997e-1,2.466405579426905e-1, 1.478980363444585e-2,3.826993559940360e-4,5.395676869878828e-6, 4.700912200921704e-8,2.733894920915608e-10,1.115830108455192e-12, 3.301093025084127e-15,7.209167098020555e-18,1.166898488777214e-20, 1.378948246502109e-23,1.124884061857506e-26,5.498556929587117e-30]
    i0q = [4.463598170691436e-1,1.702205745042606e-3, 2.792125684538934e-6,2.369902034785866e-9,8.965900179621208e-13]
    i0pp = [1.192273748120670e-1,1.947452015979746e-1, 7.629241821600588e-2,8.474903580801549e-3,2.023821945835647e-4]
    i0qq = [2.962898424533095e-1,4.866115913196384e-1, 1.938352806477617e-1,2.261671093400046e-2,6.450448095075585e-4, 1.529835782400450e-6]
    i1p = [5.000000000000000e-1,6.090824836578078e-2, 2.407288574545340e-3,4.622311145544158e-5,5.161743818147913e-7, 3.712362374847555e-9,1.833983433811517e-11,6.493125133990706e-14, 1.693074927497696e-16,3.299609473102338e-19,4.813071975603122e-22, 5.164275442089090e-25,3.846870021788629e-28,1.712948291408736e-31]
    i1q = [4.665973211630446e-1,1.677754477613006e-3, 2.583049634689725e-6,2.045930934253556e-9,7.166133240195285e-13]
    i1pp = [1.286515211317124e-1,1.930915272916783e-1, 6.965689298161343e-2,7.345978783504595e-3,1.963602129240502e-4]
    i1qq = [3.309385098860755e-1,4.878218424097628e-1, 1.663088501568696e-1,1.473541892809522e-2,1.964131438571051e-4, -1.034524660214173e-6]

    function i0(x::Float64)
        ax = abs(x)
        if ax < 15.0
            y = x * x
            return poly(i0p, 13, y) / poly(i0q, 4, 225.0 - y)
        else
            z = 1.0 - 15.0 / ax
            return exp(ax) * poly(i0pp, 4, z) / (poly(i0qq, 5, z) * sqrt(ax))
        end
    end

    function i1(x::Float64)
        ax = abs(x)
        if ax < 15.0
            y = x * x
            return x * poly(i1p, 13, y) / poly(i1q, 4, 225.0 - y)
        else
            z = 1.0 - 15.0 / ax
            ans = exp(ax) * poly(i1pp, 4, z) / (poly(i1qq, 5, z) * sqrt(ax))
            return x > 0.0 ? ans : -ans
        end
    end

    function bessel_I_n_asympt_xlarge(x::Number)
        return exp(x) / sqrt(2.0 * pi * x)
    end

    function bessel_I_n_asympt_xsmall(x::Number, n::Int)
        return ((x / 2.0)^n) / factorial(n)
    end
    




    # Explicitly for the ones we need already
    function exp_times_i0(x::Number, exp_pref::Number)
        #ax = abs(x)  # TODO check
        #if ax < 1E-1
        #    return 1
        if abs(x) < 10.0
            y = x * x
            return exp(exp_pref) * poly(i0p, 13, y) / poly(i0q, 4, 225.0 - y)
        else
            z = 1.0 - 15.0 / x
            return exp(x + exp_pref) * poly(i0pp, 4, z) / (poly(i0qq, 5, z) * sqrt(x))
        end
    end

    function exp_times_i1(x::Number, exp_pref::Number)
        #ax = abs(x)  # TODO check
        #if ax < 1E-1
        #    return (x / 2.0)
        if abs(x) < 10.0
            y = x * x
            return exp(exp_pref) * x * poly(i1p, 13, y) / poly(i1q, 4, 225.0 - y)
        else
            z = 1.0 - 15.0 / x
            return exp(x + exp_pref) * poly(i1pp, 4, z) / (poly(i1qq, 5, z) * sqrt(x))
        end
    end

    function i2(x::Number, i0_at_x::Number, i1_at_x::Number)
        if x==0
            return 0
        end

        return -2.0 / x * i1_at_x + i0_at_x
    end
end