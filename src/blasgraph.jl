mutable struct BLASGraph{T} <: AbstractBLASGraph{T}
    A::GrB_Matrix{T}
    nv::UInt64
    ne::UInt64
end
Base.show(io::IO, g::BLASGraph) = print("{", nv(g), ", ", ne(g) , "} undirected graph")

function BLASGraph(A::GrB_Matrix{T}) where T
    nrows, ncols = size(A)

    nrows != ncols && error("Matrix must be square")

    # remove all zero weight edges
    dropzeros!(A)

    # check if matrix is symmetric
    A_TRAN = A'
    if A != A_TRAN
        OK(GrB_free(A_TRAN))
        error("Matrix must be symmetric")
    end
    OK( GrB_free(A_TRAN) )

    return BLASGraph{T}(A, nrows, nnz(A)/2)
end

function BLASGraph(lg::SimpleGraph)
    A = GrB_Matrix(Int64, nv(lg), nv(lg))
    g = BLASGraph(A)
    for u in vertices(lg)
        for v in outneighbors(lg, u)
            A[OneBasedIndex(u), OneBasedIndex(v)] = 1
        end
    end
    g.ne = ne(lg)
    return g
end

is_directed(::BLASGraph) = false
is_directed(::Type{<:BLASGraph}) = false

function add_edge!(g::BLASGraph{T}, s::Union{Int64, UInt64}, d::Union{Int64, UInt64}, weight::T = one(T)) where T
    (has_edge(g, s, d) || weight == 0) && return false    # zero weight edges not allowed
    M = g.A
    u = OneBasedIndex(s)
    v = OneBasedIndex(d)
    M[u, v] = weight
    M[v, u] = weight
    g.ne += 1
    return true
end

function rem_edge!(g::BLASGraph{T}, s::Union{Int64, UInt64}, d::Union{Int64, UInt64}) where T
    has_edge(g, s, d) || return false
    M = g.A
    u = OneBasedIndex(s)
    v = OneBasedIndex(d)
    M[u, v] = T(0)
    M[v, u] = T(0)
    w = GrB_Vector(UInt64, nv(g))
    desc1 = GrB_Descriptor(Dict(GrB_OUTP => GrB_REPLACE))
    desc2 = GrB_Descriptor(Dict(GrB_OUTP => GrB_REPLACE, GrB_INP0 => GrB_TRAN))
    OK( GrB_Col_extract(w, GrB_NULL, GrB_NULL, M, GrB_ALL, 0, u, GrB_NULL) )
    OK( GrB_Col_assign(M, w, GrB_NULL, w, GrB_ALL, 0, u, desc1) )
    OK( GrB_Col_extract(w, GrB_NULL, GrB_NULL, M, GrB_ALL, 0, u, desc2) )
    OK( GrB_Row_assign(M, w, GrB_NULL, w, u, GrB_ALL, 0, desc1) )
    g.ne -= 1
    OK( GrB_free(w) )
    OK( GrB_free(desc1) )
    OK( GrB_free(desc2) )
    return true
end

function edges(g::BLASGraph)
    I, J, _ = findnz(g.A)
    e = Vector{SimpleEdge{UInt64}}(undef, ne(g))
    j::Int64 = 1
    for i = 1:length(I)
        u = I[i].x + 1
        v = J[i].x + 1
        if u <= v
            e[j] = SimpleEdge(u, v)
            j += 1
        end
    end
    return e
end

outdegree(g::BLASGraph, v::UInt64) = indegree(g, v)
outneighbors(g::BLASGraph, v::Union{Int64, UInt64}) = inneighbors(g, v)
