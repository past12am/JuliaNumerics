module NewtonRaphson
    
    function find_roots(f_df_tuple::Function, x0::Real, acc::Real)
        dx::Float64 = 1E5

        root = x0
        while dx > acc
            f_val, df_val = f_df_tuple(root)

            dx = f_val / df_val
            root -= dx
        end

        return root
    end
end