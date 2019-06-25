indexing(n) = [(i,j) for i in 1:n for j in 1:n if i≠j]

function generating_set(G::AutGroup{N}, n=N) where N

    rmuls = [Groups.rmul_autsymbol(i,j) for (i,j) in indexing(n)]
    lmuls = [Groups.lmul_autsymbol(i,j) for (i,j) in indexing(n)]
    gen_set = G.([rmuls; lmuls])

    return [gen_set; inv.(gen_set)]
end

function E(M::MatSpace, i::Integer, j::Integer, val=1)
    @assert i ≠ j
    @assert 1 ≤ i ≤ nrows(M)
    @assert 1 ≤ j ≤ ncols(M)
    m = one(M)
    m[i,j] = val
    return m
end

function generating_set(M::MatSpace, n=nrows(M))
    @assert nrows(M) == ncols(M)

    elts = [E(M, i,j) for (i,j) in indexing(n)]
    return elem_type(M)[elts; inv.(elts)]
end

isopposite(σ::perm, τ::perm, i=1, j=2) =
    σ[i] ≠ τ[i] && σ[i] ≠ τ[j] &&
    σ[j] ≠ τ[i] && σ[j] ≠ τ[j]

isadjacent(σ::perm, τ::perm, i=1, j=2) =
    (σ[i] == τ[i] && σ[j] ≠ τ[j]) || # first equal, second differ
    (σ[j] == τ[j] && σ[i] ≠ τ[i]) || # sedond equal, first differ
    (σ[i] == τ[j] && σ[j] ≠ τ[i]) || # first σ equal to second τ
    (σ[j] == τ[i] && σ[i] ≠ τ[j])    # second σ equal to first τ

function Sq(RG::GroupRing, N::Integer)
    S₂ = generating_set(RG.group, 2)
    ℤ = Int64
    Δ₂ = length(S₂)*one(RG, ℤ) - RG(S₂, ℤ);

    Alt_N = [g for g in PermutationGroup(N) if parity(g) == 0]

    sq = RG()
    for σ in Alt_N
        GroupRings.addeq!(sq, *(σ(Δ₂), σ(Δ₂), false))
    end
    return RG(sq.coeffs.÷factorial(N-2))
end

function Adj(RG::GroupRing, N::Integer)
    S₂ = generating_set(RG.group, 2)
    ℤ = Int64
    Δ₂ = length(S₂)*one(RG, ℤ) - RG(S₂, ℤ);

    Alt_N = [g for g in PermutationGroup(N) if parity(g) == 0]
    Δ₂s = Dict(σ=>σ(Δ₂) for σ in Alt_N)
    adj = RG()

    for σ in Alt_N
        for τ in Alt_N
            if isadjacent(σ, τ)
                GroupRings.addeq!(adj, *(Δ₂s[σ], Δ₂s[τ], false))
            end
        end
    end
    return RG(adj.coeffs.÷factorial(N-2)^2)
end

function Op(RG::GroupRing, N::Integer)
    if N < 4
        return RG()
    end
    S₂ = generating_set(RG.group, 2)
    ℤ = Int64
    Δ₂ = length(S₂)*one(RG, ℤ) - RG(S₂, ℤ);

    Alt_N = [g for g in PermutationGroup(N) if parity(g) == 0]
    Δ₂s = Dict(σ=>σ(Δ₂) for σ in Alt_N)
    op = RG()

    for σ in Alt_N
        for τ in Alt_N
            if isopposite(σ, τ)
                GroupRings.addeq!(op, *(Δ₂s[σ], Δ₂s[τ], false))
            end
        end
    end
    return RG(op.coeffs.÷factorial(N-2)^2)
end

function Ygrek(RG::GroupRing, N)
    ygrek = RG()
    elt = RG(generating_set(RG.group)[1])

    G = PermutationGroup(N)
    elts = Dict(g=>g(elt) for g in G)

    for g in G
        for h in G
            if isadjacent(g,h)
                a = elts[g] - elts[h]
                GroupRings.addeq!(ygrek, *(a, a, false))
                a = star(elts[g]) - elts[h]
                GroupRings.addeq!(ygrek, *(a, a, false))
                a = star(elts[g]) - star(elts[h])
                GroupRings.addeq!(ygrek, *(a, a, false))
                a = elts[g] - star(elts[h])
                GroupRings.addeq!(ygrek, *(a, a, false))
            end
        end
    end
    return RG(ygrek.coeffs.÷(2factorial(N-2)^2))
end


for Op in [:Sq, :Adj, :Op, :Ygrek]
    @eval begin
        $Op(RG::GroupRing{AutGroup{N}}) where N = $Op(RG, N)
        $Op(RG::GroupRing{<:MatSpace}) = $Op(RG, nrows(RG.group))
    end
end

function SqAdjOp(RG::GroupRing, N::Integer)
    S₂ = generating_set(RG.group, 2)
    ℤ = Int64
    Δ₂ = length(S₂)*one(RG, ℤ) - RG(S₂, ℤ);

    Alt_N = [σ for σ in PermutationGroup(N) if parity(σ) == 0]
    sq, adj, op = RG(), RG(), RG()

    Δ₂s = Dict(σ=>σ(Δ₂) for σ in Alt_N)

    for σ in Alt_N
        GroupRings.addeq!(sq, *(Δ₂s[σ], Δ₂s[σ], false))
        for τ in Alt_N
            if isopposite(σ, τ)
                GroupRings.addeq!(op, *(Δ₂s[σ], Δ₂s[τ], false))
            elseif isadjacent(σ, τ)
                GroupRings.addeq!(adj, *(Δ₂s[σ], Δ₂s[τ], false))
            end
        end
    end

    k = factorial(N-2)
    return RG(sq.coeffs.÷k), RG(adj.coeffs.÷k^2), RG(op.coeffs.÷k^2)
end
