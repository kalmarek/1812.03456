using SparseArrays

using AbstractAlgebra
using Groups
using GroupRings
using PropertyT

using IntervalArithmetic

using SCS
using JLD

include(joinpath("src", "sqadjop.jl"))

invalid_use_message = """You need to call this script in the parent folder of oSAutF5_r2 folder.
Provide also the two (numerical) parameters: `-k` and `-lambda`"""

if !(iseven(length(ARGS)))
    throw(invalid_use_message)
end

K = LAMBDA = nothing 

for i in 1:2:length(ARGS)
    if ARGS[i] == "-k"
        K = parse(Float64, ARGS[i+1])
    elseif ARGS[i] == "-lambda"
        LAMBDA = parse(Float64, ARGS[i+1])
    end
end

if K == nothing || LAMBDA == nothing
    throw(invalid_use_message)
end

@info("Running checks for Adj₅ + $K·Op₅ - $LAMBDA·Δ₅") 

N = 5
G = AutGroup(FreeGroup(N), special=true)
S = generating_set(G)

const prefix = "oSAutF$(N)_r2"
isdir(prefix) || mkpath(prefix)

const DELTA_FILE = joinpath(prefix,"delta.jld")
const SQADJOP_FILE = joinpath(prefix, "SqAdjOp_coeffs.jld")
const ORBITDATA_FILE = joinpath(prefix, "OrbitData.jld")

const fullpath = joinpath(prefix, string(LAMBDA))
isdir(fullpath) || mkpath(fullpath)
const SOLUTION_FILE = joinpath(fullpath, "solution.jld")

@info("Looking for delta.jld, SqAdjOp_coeffs.jld and OrbitData.jld in $prefix")

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
    info("Computing Laplacian")
    Δ = PropertyT.Laplacian(S, 2)
    PropertyT.saveGRElem(DELTA_FILE, Δ)
    RG = parent(Δ)

    info("Computing Sq, Adj, Op")
    @time sq, adj, op = Sq(RG), Adj(RG), Op(RG)
    
    save(SQADJOP_FILE, "Sq", sq.coeffs, "Adj", adj.coeffs, "Op", op.coeffs)

    info("Compute OrbitData")
    if !isfile(ORBITDATA_FILE)
        orbit_data = PropertyT.OrbitData(RG, sett.autS)
        save(ORBITDATA_FILE, "OrbitData", orbit_data)
    else
        orbit_data = load(ORBITDATA_FILE, "OrbitData")
    end
end;

orbit_data = PropertyT.decimate(orbit_data);

elt = adj + K*op;

info("Looking for solution.jld in $fullpath")

if !isfile(SOLUTION_FILE)
    info("solution.jld not found, attempting to recreate one.")
    
    SDP_problem, varλ, varP = PropertyT.SOS_problem(elt, Δ, orbit_data; upper_bound=LAMBDA)
    
    begin
        scs_solver = SCS.SCSSolver(linear_solver=SCS.Direct,
            eps=1e-12,
            max_iters=200_000,
            alpha=1.5,
            acceleration_lookback=1)

        JuMP.setsolver(SDP_problem, scs_solver)
    end
    
    λ = Ps = ws = nothing
    const WARMSTART_FILE = joinpath(fullpath, "warmstart.jld")
    if isfile(WARMSTART_FILE)
        ws = load(WARMSTART_FILE, "warmstart")
    end

    i = 0
    # for i in 1:6
    status= :Unknown
    while status !=:Optimal
        i += 1
        SOLVERLOG_FILE = joinpath(fullpath, "solver_$(now()).log")
        status, (λ, Ps, ws) = PropertyT.solve(SOLVERLOG_FILE, SDP_problem, varλ, varP, ws);
        precision = abs(λ - LAMBDA)
        println("i = $i, \t precision = $precision")
        
        if all((!isnan).(ws[1]))
            save(WARMSTART_FILE, "warmstart", ws, "λ", λ, "Ps", Ps)
            save(WARMSTART_FILE[1:end-4]*"_$(now()).jld", "warmstart", ws, "λ", λ, "Ps", Ps)
        else
            warn("No valid solution was saved!")
        end
    end

    info("Reconstructing P...")
    @time P = PropertyT.reconstruct(Ps, orbit_data);
    info("Computing Q = √P")
    @time const Q = real(sqrtm(P));
    
    save(SOLUTION_FILE, "λ", λ, "Q", Q)
end

info("Checking the sum of squares solution for 36(Adj₅ + $K·Op₅) - $LAMBDA·Δ₅") 
Q, λ = load(SOLUTION_FILE, "Q", "λ")

function SOS_residual(eoi::GroupRingElem, Q::Matrix)
    RG = parent(eoi)
    @time sos = PropertyT.compute_SOS(RG, Q);
    return eoi - sos
end

info("In floating point arithmetic:")
EOI = elt - λ*Δ
b = SOS_residual(EOI, Q)
@show norm(b, 1);

info("In interval arithmetic:")
EOI_int = elt - @interval(λ)*Δ;
Q_int = PropertyT.augIdproj(Q);
@assert all([zero(eltype(Q)) in sum(view(Q_int, :, i)) for i in 1:size(Q_int, 2)])
b_int = SOS_residual(EOI_int, Q_int)
@show norm(b_int, 1);

λ_cert = @interval(λ) - 2^2*norm(b_int,1)
info("λ is certified to be > ", λ_cert.lo)

info("i.e Adj₅ + $K·Op₅ - $((λ_cert/36).lo)·Δ₅ ∈ Σ²₂ ISAut(F₅)")
