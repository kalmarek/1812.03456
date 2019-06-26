function parse_args(args)
    invalid_use_message = """You need to call this script in the parent folder of oSAutF5_r2 folder.
Provide also the two (numerical) parameters: `--k` and `--lambda`"""

    iseven(length(args)) || throw(invalid_use_message)

    k, λ = nothing, nothing

    for i in 1:2:length(ARGS)
        if ARGS[i] == "--k"
            k = try
                parse(Float64, ARGS[i+1])
            catch
                throw(invalid_use_message)
            end
        elseif ARGS[i] == "--lambda"
            λ = try
                parse(Float64, ARGS[i+1])
            catch
                throw(invalid_use_message)
            end
        end
    end
    @show k, λ
    
    if isnothing(k) || isnothing(λ)
        throw(invalid_use_message)
    end
    return k, λ
end 
