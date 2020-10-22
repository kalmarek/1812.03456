using ArgParse

function parse_commandline()
    s = ArgParseSettings(
    description="This program computes and/or verifies an approximate sum of squares decomposition for element
    `Adj_n + k·Op_n - l·Δ_n ~ Σξ_i*ξ_i ∈ Σ²₂ ISAut(F_n)`.
    Parameters `-n`, `-k` and `-l` must be provided when executing the script. \n",
    preformatted_description=true)
    @add_arg_table! s begin
        "-n"
            help = "Perform computations in automorphism group of the free group on `n` generators."
            arg_type = Int
            required = true
        "-k"
            help = "Compute sum of squares for `Adj_n + k·Op_n ∈ ISAut(F_n)`."
            arg_type = Float64
            required = true
        "--lambda", "-l"
            help = "upper bound on the amount of `Δ_n` contained in the element."
            arg_type = Float64
            required = true
        "--tol"
            help = "solvers numerical tolerance"
            required = false
            default = 1e-12
    end

    return parse_args(s)
end
