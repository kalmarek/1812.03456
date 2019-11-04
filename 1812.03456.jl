using SparseArrays
using LinearAlgebra
using Dates

using AbstractAlgebra
using Groups
using GroupRings
using PropertyT
using SCS

using PropertyT.JuMP
using PropertyT.IntervalArithmetic

using PropertyT.JLD

include(joinpath("src", "argparse.jl"))
N, K, LAMBDA, scseps = parse_args(ARGS)

@info "Running checks for Adj_$N + $K·Op_$N - $LAMBDA·Δ_$N"

G = AutGroup(FreeGroup(N), special=true)
S = PropertyT.generating_set(G)

const prefix = "oSAutF$(N)_r2"
isdir(prefix) || mkpath(prefix)

const DELTA_FILE = joinpath(prefix,"delta.jld")
const SQADJOP_FILE = joinpath(prefix, "SqAdjOp_coeffs.jld")
const ORBITDATA_FILE = joinpath(prefix, "OrbitData.jld")

const fullpath = joinpath(prefix, "$(LAMBDA)_K=$K")
isdir(fullpath) || mkpath(fullpath)
const WARMSTART_FILE = joinpath(fullpath, "warmstart.jld")
const SOLUTION_FILE = joinpath(fullpath, "solution.jld")

@info "Looking for delta.jld, SqAdjOp_coeffs.jld and OrbitData.jld in $prefix"

if isfile(DELTA_FILE) && isfile(SQADJOP_FILE) && isfile(ORBITDATA_FILE)
    # cached
    Δ = PropertyT.loadGRElem(DELTA_FILE, G)
    RG = parent(Δ)
    orbit_data = load(ORBITDATA_FILE, "OrbitData")
    sq_c, adj_c, op_c = load(SQADJOP_FILE, "Sq", "Adj", "Op")
    sq = GroupRingElem(sq_c, RG)
    adj = GroupRingElem(adj_c, RG)
    op = GroupRingElem(op_c, RG);
else
    @info "Computing Laplacian"
    Δ = PropertyT.Laplacian(S, 2)
    PropertyT.saveGRElem(DELTA_FILE, Δ)
    RG = parent(Δ)

    @info "Computing Sq, Adj, Op"
    @time sq, adj, op = PropertyT.SqAdjOp(RG, N)

    save(SQADJOP_FILE, "Sq", sq.coeffs, "Adj", adj.coeffs, "Op", op.coeffs)

    @info "Compute OrbitData"
    if !isfile(ORBITDATA_FILE)
        orbit_data = PropertyT.OrbitData(RG, WreathProduct(PermGroup(2), PermGroup(N)))
        save(ORBITDATA_FILE, "OrbitData", orbit_data)
    else
        orbit_data = load(ORBITDATA_FILE, "OrbitData")
    end
end;

orbit_data = PropertyT.decimate(orbit_data);

elt = adj + K*op;
ELT_STRING = "Adj_$(N)+$(K)·Op_$(N)"

@info "Looking for solution.jld in $fullpath"

interrupted = false

if !isfile(SOLUTION_FILE)
    @info "$SOLUTION_FILE not found, attempting to recreate one."

    SDP_problem, varP = PropertyT.SOS_problem(elt, Δ, orbit_data; upper_bound=LAMBDA)

    with_SCS = JuMP.with_optimizer(SCS.Optimizer, linear_solver=SCS.Direct,
                             max_iters=500_000,
                             eps=scseps,
                             alpha=1.5,
                             acceleration_lookback=1,
                             warm_start=true)

    let
        λ = Ps = ws = nothing
        if isfile(WARMSTART_FILE)
            ws = load(WARMSTART_FILE, "warmstart")
        end

        interrupted = false
        status = nothing
        while status != JuMP.MOI.OPTIMAL
            SOLVERLOG_FILE = joinpath(fullpath, "$(ELT_STRING)_solver_$(now()).log")
            @info "Recording solvers progress in" SOLVERLOG_FILE
            @time status, ws = PropertyT.solve(SOLVERLOG_FILE, SDP_problem, with_SCS, ws);
            λ = value(SDP_problem[:λ])
            Ps = [value.(P) for P in varP]

            if all((!isnan).(ws[1])) # solution looks valid
                save(WARMSTART_FILE,
                    "warmstart", (ws.primal, ws.dual, ws.slack), "Ps", Ps, "λ", λ)
                save(WARMSTART_FILE[1:end-4]*"$(now())"*".jld",
                    "warmstart", (ws.primal, ws.dual, ws.slack), "Ps", Ps, "λ", λ)
            else
                @warn "No valid solution was saved!"
                interrupted = true
                break
            end
        end

        if interrupted
            @info "Reading the last saved warmstart for λ and Ps"
            λ, Ps = load(WARMSTART_FILE, "λ", "Ps")
        end

        @info "Reconstructing Q..."
        Qs = real.(sqrt.(Ps))
        @time Q = PropertyT.reconstruct(Qs, orbit_data);
        save(SOLUTION_FILE, "λ", λ, "Q", Q)
    end
end

@info "Checking the sum of squares solution for Adj_N + $K·Op_N - $LAMBDA·Δ_N"
Q, λ = load(SOLUTION_FILE, "Q", "λ")

function SOS_residual(eoi::GroupRingElem, Q::Matrix)
    RG = parent(eoi)
    @time sos = PropertyT.compute_SOS(RG, Q);
    return eoi - sos
end

let EOI = elt - λ*Δ
    residual = SOS_residual(EOI, Q)
    @info "In floating point arithmetic the ℓ₁-norm of the residual" norm(residual, 1);
end

let EOI_int = elt - @interval(λ)*Δ;
    Q_int, check_columns_in_augmentation_ideal = PropertyT.augIdproj(Interval, Q);
    @assert check_columns_in_augmentation_ideal

    residual_int = SOS_residual(EOI_int, Q_int);
    @info "In interval arithmetic the ℓ₁-norm of the residual" norm(residual_int, 1);

    λ_cert = @interval(λ) - 2^2*norm(residual_int,1)
    @info "λ is certified to be > " λ_cert.lo
    @info "i.e Adj_$N + $K·Op_$N - ($(λ_cert.lo))·Δ_$N ∈ Σ²₂ ISAut(F_$N)"
end
