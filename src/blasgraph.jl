mutable struct BLASGraph{T} <: AbstractBLASGraph{T}
    A::GrB_Matrix{T}
    nv::UInt64
    ne::UInt64
end
Base.show(io::IO, g::BLASGraph) = print("{", nv(g), ", ", ne(g) , "} undirected graph")

BLASGraph{T}(A::GrB_Matrix{T}) where T = BLASGraph(A)

function BLASGraph(A::GrB_Matrix{T}) where T
    nrows, ncols = size(A)

    nrows != ncols && error("Matrix must be square")

    # remove all zero weight edges
    SuiteSparseGraphBLAS.dropzeros!(A)

    # check if matrix is symmetric
    A_TRAN = A'
    if A != A_TRAN
        OK(GrB_free(A_TRAN))
        error("Matrix must be symmetric")
    end
    OK( GrB_free(A_TRAN) )

    return BLASGraph{T}(A, nrows, SuiteSparseGraphBLAS.nnz(A)/2)
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

function BLASGraph(edge_list::Array{SimpleWeightedEdge{Int64, T}, 1}) where T
    nvg = 0
    for e in edge_list
        nvg = max(nvg, src(e), dst(e))
    end
    g = BLASGraph(GrB_Matrix(T, nvg, nvg))
    for e in edge_list
        add_edge!(g, e)
    end
    return g
end

BLASGraph(n::Union{Int64, UInt64}) = BLASGraph(GrB_Matrix(Float64, n, n))
BLASGraph{T}(n::Union{Int64, UInt64}) where T = BLASGraph(GrB_Matrix(T, n, n))

function BLASGraph(adjmx::SparseMatrixCSC)
    I, J, X = SparseArrays.findnz(adjmx)
    return BLASGraph(GrB_Matrix(OneBasedIndex.(I), OneBasedIndex.(J), X, nrows = size(adjmx, 1), ncols = size(adjmx, 2)))
end

BLASGraph(m::AbstractMatrix) = BLASGraph(sparse(m))

is_directed(::BLASGraph) = false
is_directed(::Type{<:BLASGraph}) = false

function add_edge!(g::BLASGraph{T}, s::Integer, d::Integer, weight::T = one(T)) where T
    (!has_vertex(g, s) || !has_vertex(g, d) || weight == 0) && return false    # zero weight edges not allowed
    edge_was_present = has_edge(g, s, d)
    M = g.A
    u = OneBasedIndex(s)
    v = OneBasedIndex(d)
    M[u, v] = weight
    M[v, u] = weight
    g.ne += 1
    return true
end

function rem_edge!(g::BLASGraph{T}, s::Integer, d::Integer) where T
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

function edges(g::BLASGraph{T}) where T
    I, J, X = SuiteSparseGraphBLAS.findnz(g.A)
    e = Vector{SimpleWeightedEdge{UInt64, T}}(undef, ne(g))
    j = 1
    for i = 1:length(I)
        u = I[i].x + 1
        v = J[i].x + 1
        if u <= v
            w = X[i]
            e[j] = SimpleWeightedEdge(u, v, w)
            j += 1
        end
    end
    return e
end

outdegree(g::BLASGraph, v::Integer) = indegree(g, v)
outneighbors(g::BLASGraph, v::Integer) = inneighbors(g, v)
