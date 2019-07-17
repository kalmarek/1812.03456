function parse_args(args)
    invalid_use_message = """You need to call this script in the parent folder of oSAutF5_r2 folder.
Provide also the two (numerical) parameters: `--k` (`-k`) and `--lambda` (`-l`).
Optional prameters:
\t `--n` (`-n`) [5]     : check/reconstruct solution for SAut(F_n)
\t `--eps`	[1e-12] : SCS solver precision
"""

    iseven(length(args)) || throw(invalid_use_message)

    n, k, λ = 5, nothing, nothing, 1e-12

    for i in 1:2:length(ARGS)
        arg = ARGS[i]
        next_arg = ARGS[i+1]
        if arg == "--k" || arg == "-k"
            k = try
                parse(Float64, next_arg)
            catch
                throw(invalid_use_message)
            end
        elseif arg == "--lambda" || arg == "-l"
            λ = try
                parse(Float64, next_arg)
            catch
                throw(invalid_use_message)
            end
        elseif arg == "--n" || arg == "-n"
            n = try
                parse(Int, next_arg)
            catch
                throw(invalid_use_message)
            end
        elseif arg == "--eps"
            eps = try
                parse(Float64, next_arg)
            catch
                throw(invalid_use_message)
            end
        end
    end
    
    if isnothing(k) || isnothing(λ)
        throw(invalid_use_message)
    end
    return n, k, λ, eps
end 
