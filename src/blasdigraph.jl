mutable struct BLASDiGraph{T} <: AbstractBLASGraph{T}
    A::GrB_Matrix{T}
    nv::UInt64
    ne::UInt64
end
Base.show(io::IO, g::BLASDiGraph) = print("{", nv(g), ", ", ne(g), "} directed graph")

BLASDiGraph{T}(A::GrB_Matrix{T}) where T = BLASDiGraph(A)

function BLASDiGraph(A::GrB_Matrix{T}) where T
    nrows, ncols = size(A)

    nrows != ncols && error("Matrix must be square")

    # remove all zero weight edges
    SuiteSparseGraphBLAS.dropzeros!(A)

    return BLASDiGraph{T}(A, nrows, SuiteSparseGraphBLAS.nnz(A))
end

function BLASDiGraph(lg::SimpleDiGraph)
    A = GrB_Matrix(Int64, nv(lg), nv(lg))
    g = BLASDiGraph(A)
    for u in vertices(lg)
        for v in outneighbors(lg, u)
            A[OneBasedIndex(u), OneBasedIndex(v)] = 1
        end
    end
    g.ne = ne(lg)
    return g
end

BLASDiGraph(n::Union{Int64, UInt64}) = BLASDiGraph(GrB_Matrix(Float64, n, n))
BLASDiGraph{T}(n::Union{Int64, UInt64}) where T = BLASDiGraph(GrB_Matrix(T, n, n))

function BLASDiGraph(adjmx::SparseMatrixCSC)
    I, J, X = SparseArrays.findnz(adjmx)
    return BLASDiGraph(GrB_Matrix(OneBasedIndex.(I), OneBasedIndex.(J), X, nrows = size(adjmx, 1), ncols = size(adjmx, 2)))
end

BLASDiGraph(m::AbstractMatrix) = BLASDiGraph(sparse(m))

function outdegree(g::BLASDiGraph, v::Integer)
    M = g.A
    row_v = GrB_Vector(Bool, nv(g))
    inp0_tran_desc = GrB_Descriptor(Dict(GrB_INP0 => GrB_TRAN))
    OK( GrB_Col_extract(row_v, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), OneBasedIndex(v), inp0_tran_desc) )
    n = SuiteSparseGraphBLAS.nnz(row_v)
    OK( GrB_free(row_v) )
    OK( GrB_free(inp0_tran_desc) )
    return n
end

is_directed(::BLASDiGraph) = true
is_directed(::Type{<:BLASDiGraph}) = true

function add_edge!(g::BLASDiGraph{T}, s::Integer, d::Integer, weight::T = one(T)) where T
    (!has_vertex(g, s) || !has_vertex(g, d) || weight == 0) && return false    # zero weight edges not allowed
    M = g.A
    M[OneBasedIndex(s), OneBasedIndex(d)] = weight
    g.ne += 1
    return true
end

function rem_edge!(g::BLASDiGraph{T}, s::Integer, d::Integer) where T
    has_edge(g, s, d) || return false
    M = g.A
    u = OneBasedIndex(s)
    v = OneBasedIndex(d)
    M[u, v] = T(0)
    w = GrB_Vector(UInt64, nv(g))
    desc1 = GrB_Descriptor(Dict(GrB_OUTP => GrB_REPLACE))
    OK( GrB_Col_extract(w, GrB_NULL, GrB_NULL, M, GrB_ALL, 0, u, GrB_NULL) )
    OK( GrB_Col_assign(M, w, GrB_NULL, w, GrB_ALL, 0, u, desc1) )
    g.ne -= 1
    OK( GrB_free(w) )
    return true
end

function edges(g::BLASDiGraph{T}) where T
    I, J, X = SuiteSparseGraphBLAS.findnz(g.A)
    e = Vector{SimpleWeightedEdge{UInt64, T}}(undef, ne(g))
    for i = 1:length(I)
        u = I[i].x + 1
        v = J[i].x + 1
        w = X[i]
        e[i] = SimpleWeightedEdge(u, v, w)
    end
    return e
end

function outneighbors(g::BLASDiGraph, v::Integer)
    M = g.A
    x = GrB_Vector(Bool, nv(g))
    inp0_tran_desc = GrB_Descriptor(Dict(GrB_INP0 => GrB_TRAN))
    OK( GrB_Col_extract(x, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), OneBasedIndex(v), inp0_tran_desc) )
    I, X = SuiteSparseGraphBLAS.findnz(x)
    nbrs = Vector{UInt64}(undef, length(I))
    for i = 1:length(I)
        nbrs[i] = I[i].x + 1
    end
    OK( GrB_free(x) )
    OK( GrB_free(inp0_tran_desc) )
    return nbrs
end
