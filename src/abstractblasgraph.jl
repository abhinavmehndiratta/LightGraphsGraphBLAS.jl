abstract type AbstractBLASGraph{T} <: AbstractGraph{UInt64} end

ne(g::AbstractBLASGraph) = g.ne
nv(g::AbstractBLASGraph) = g.nv

function indegree(g::AbstractBLASGraph, v::Union{Int64, UInt64})
    M = g.A
    col_v = GrB_Vector(Bool, nv(g))
    OK( GrB_Col_extract(col_v, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), OneBasedIndex(v), GrB_NULL) )
    n = nnz(col_v)
    OK( GrB_free(col_v) )
    return n
end

vertices(g::AbstractBLASGraph) = 1:nv(g)

has_vertex(g::AbstractBLASGraph, v::Union{Int64, UInt64}) = (v >= 1 && v <= nv(g)) ? true : false

function has_edge(g::AbstractBLASGraph, s::Union{Int64, UInt64}, d::Union{Int64, UInt64})
    M = g.A
    (typeof(GrB_Matrix_extractElement(M, OneBasedIndex(s), OneBasedIndex(d))) != GrB_Info) && return true
    return false
end

function inneighbors(g::AbstractBLASGraph, v::Union{Int64, UInt64})
    M = g.A
    x = GrB_Vector(Bool, nv(g))
    OK( GrB_Col_extract(x, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), OneBasedIndex(v), GrB_NULL) )
    I, X = findnz(x)
    OK( GrB_free(x) )
    nbrs = Vector{UInt64}(undef, length(I))
    for i = 1:length(I)
        nbrs[i] = I[i].x + 1
    end
    return nbrs
end

eltype(::AbstractBLASGraph) = UInt64

edgetype(::AbstractBLASGraph) = SimpleEdge{UInt64}

weights(g::AbstractBLASGraph) = BLASGraphWeights(g.A)
