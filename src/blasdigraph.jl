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

function BLASDiGraph(edge_list::Array{SimpleWeightedEdge{Int64, T}, 1}) where T
    nvg = 0
    for e in edge_list
        nvg = max(nvg, src(e), dst(e))
    end
    g = BLASDiGraph(GrB_Matrix(T, nvg, nvg))
    for e in edge_list
        add_edge!(g, e)
    end
    return g
end

BLASDiGraph(n::Union{Int64, UInt64}) = BLASDiGraph(GrB_Matrix(Float64, n, n))
BLASDiGraph{T}(n::Union{Int64, UInt64}) where T = BLASDiGraph(GrB_Matrix(T, n, n))

function BLASDiGraph(adjmx::AbstractMatrix{T}) where T
    dima, dimb = size(adjmx)

    dima != dimb && error("Matrix must be square")

    g = BLASDiGraph{T}(dima)
    A = g.A

    for i in findall(adjmx .!= 0)
        r = i[1]
        c = i[2]
        w = adjmx[r, c]
        A[OneBasedIndex(r), OneBasedIndex(c)] = w
        g.ne += 1
    end

    return g
end

function BLASDiGraph(adjmx::SparseMatrixCSC{T}) where T
    dima, dimb = size(adjmx)
    dima != dimb && error("Matrix must be square")

    g = BLASDiGraph{T}(dima)
    A = g.A
    maxc = length(adjmx.colptr)

    for c = 1:(maxc - 1)
        for rind = adjmx.colptr[c]:(adjmx.colptr[c + 1] - 1)
            w = adjmx.nzval[rind]
            if w != 0
                r = adjmx.rowval[rind]
                A[OneBasedIndex(r), OneBasedIndex(c)] = w
                g.ne += 1
            end
        end
    end

    return g
end

function outdegree(g::BLASDiGraph, v::Integer)
    M = g.A
    row_v = GrB_Vector(Bool, nv(g))
    inp0_tran_desc = GrB_Descriptor(Dict(GrB_INP0 => GrB_TRAN))
    OK( GrB_Col_extract(row_v, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), OneBasedIndex(v), inp0_tran_desc) )
    n = SuiteSparseGraphBLAS.nnz(row_v)
    OK( GrB_free(row_v) )
    OK( GrB_free(inp0_tran_desc) )
    return Int64(n)
end

is_directed(::BLASDiGraph) = true
is_directed(::Type{<:BLASDiGraph}) = true

function add_edge!(g::BLASDiGraph{T}, s::Integer, d::Integer, weight::T = one(T)) where T
    (!has_vertex(g, s) || !has_vertex(g, d) || weight == 0) && return false    # zero weight edges not allowed
    edge_was_present = has_edge(g, s, d)
    M = g.A
    M[OneBasedIndex(s), OneBasedIndex(d)] = weight
    if !edge_was_present
        g.ne += 1
    end
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
    desc2 = GrB_Descriptor(Dict(GrB_OUTP => GrB_REPLACE, GrB_INP0 => GrB_TRAN))
    OK( GrB_Col_extract(w, GrB_NULL, GrB_NULL, M, GrB_ALL, 0, u, desc2) )
    OK( GrB_Row_assign(M, w, GrB_NULL, w, u, GrB_ALL, 0, desc1) )
    g.ne -= 1
    OK( GrB_free(w) )
    OK( GrB_free(desc1) )
    OK( GrB_free(desc2) )
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
