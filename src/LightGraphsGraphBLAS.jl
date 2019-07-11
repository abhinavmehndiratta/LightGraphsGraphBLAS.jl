module LightGraphsGraphBLAS

using LightGraphs, GraphBLASInterface, SuiteSparseGraphBLAS

import LightGraphs:
    ne, nv, is_directed, has_vertex, vertices, has_edge, add_edge!, rem_edge!,
    indegree, outdegree, inneighbors, outneighbors, edgetype, edges,
    SimpleGraphs.SimpleEdge, weights

import Base:
    eltype, zero, size, getindex, setindex!

# check if GraphBLAS operation was successful, else throw an error
function OK(info::GrB_Info)
    info != GrB_SUCCESS && error(info)
    return true
end

abstract type AbstractBLASGraph{T} <: AbstractGraph{Int64} end

mutable struct BLASGraph{T} <: AbstractBLASGraph{T}
        A::GrB_Matrix{T}
        nv::Int64
        ne::Int64
end
Base.show(io::IO, g::BLASGraph) = print("{", nv(g), ", ", ne(g) , "} undirected graph")

mutable struct BLASDiGraph{T} <: AbstractBLASGraph{T}
        A::GrB_Matrix{T}
        nv::Int64
        ne::Int64
end
Base.show(io::IO, g::BLASDiGraph) = print("{", nv(g), ", ", ne(g), "} directed graph")

struct BLASGraphWeights{T} <: AbstractMatrix{T}
        A::GrB_Matrix{T}
end

size(w::BLASGraphWeights) = size(w.A)
getindex(w::BLASGraphWeights, r::Int64, c::Int64) = getindex(w.A, r, c)
setindex!(w::BLASGraphWeights{T}, X::T, I::Int64, J::Int64) where T = setindex!(w.A, X, I, J)

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

function BLASDiGraph(A::GrB_Matrix{T}) where T
    nrows, ncols = size(A)

    nrows != ncols && error("Matrix must be square")

    # remove all zero weight edges
    dropzeros!(A)

    return BLASDiGraph{T}(A, nrows, nnz(A))
end

function free(g::BLASGraph)
    OK( GrB_free(g.A) )
    OK( GrB_free(g.degree) )
end

function free(g::BLASDiGraph)
    OK( GrB_free(g.A) )
    OK( GrB_free(g.indegree) )
    OK( GrB_free(g.outdegree) )
end

ne(g::AbstractBLASGraph) = g.ne
nv(g::AbstractBLASGraph) = g.nv

function indegree(g::AbstractBLASGraph, v::Int64)
    M = g.A
    col_v = GrB_Vector(Bool, nv(g))
    OK( GrB_Col_extract(col_v, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), v, GrB_NULL) )
    n = nnz(col_v)
    OK( GrB_free(col_v) )
    return n
end

function outdegree(g::BLASDiGraph, v::Int64)
    M = g.A
    row_v = GrB_Vector(Bool, nv(g))
    inp0_tran_desc = GrB_Descriptor(Dict(GrB_INP0 => GrB_TRAN))
    OK( GrB_Col_extract(row_v, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), v, inp0_tran_desc) )
    n = nnz(row_v)
    OK( GrB_free(row_v) )
    OK( GrB_free(inp0_tran_desc) )
    return n
end

outdegree(g::BLASGraph) = indegree(g)

is_directed(::BLASGraph) = false
is_directed(::BLASDiGraph) = true
is_directed(::Type{<:BLASGraph}) = false
is_directed(::Type{<:BLASDiGraph}) = true

vertices(g::AbstractBLASGraph) = 0:(nv(g)-1)

has_vertex(g::AbstractBLASGraph, v::GrB_Index) = (v >= 0 && v < nv(g)) ? true : false

function add_edge!(g::BLASGraph{T}, s::GrB_Index, d::GrB_Index, weight::T = one(T)) where T
    M = g.A
    (has_edge(g, s, d) || weight == zero(T)) && return false    # zero weight edges not allowed
    try
        M[s, d] = weight
        M[d, s] = weight
        g.ne += 1
        return true
    catch
        return false
    end
end

function add_edge!(g::BLASDiGraph{T}, s::GrB_Index, d::GrB_Index, weight::T = one(T)) where T
    M = g.A
    (has_edge(g, s, d) || weight == zero(T)) && return false    # zero weight edges not allowed
    try
        M[s, d] = weight
        g.ne += 1
        return true
    catch
        return false
    end
end

function rem_edge!(g::BLASGraph{T}, s::GrB_Index, d::GrB_Index) where T
    has_edge(g, s, d) || return false
    M = g.A
    try
        M[s, d] = T(0)
        M[d, s] = T(0)
        dropzeros!(M)
        g.ne -= 1
        return true
    catch
        return false
    end
end

function rem_edge!(g::BLASDiGraph{T}, s::GrB_Index, d::GrB_Index) where T
    has_edge(g, s, d) || return false
    M = g.A
    try
        M[s, d] = T(0)
        dropzeros!(M)
        g.ne -= 1
        return true
    catch
        return false
    end
end

function has_edge(g::AbstractBLASGraph, s::GrB_Index, d::GrB_Index)
    M = g.A
    try
        M[s, d]
        return true
    catch
        return false
    end
end

function inneighbors(g::AbstractBLASGraph, v::GrB_Index)
    M = g.A
    x = GrB_Vector(Bool, nv(g))
    OK( GrB_Col_extract(x, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), v, GrB_NULL) )
    nbrs, _ = findnz(x)
    OK( GrB_free(x) )
    return nbrs
end

function outneighbors(g::AbstractBLASGraph, v::GrB_Index)
    M = g.A
    x = GrB_Vector(Bool, nv(g))
    inp0_tran_desc = GrB_Descriptor(Dict(GrB_INP0 => GrB_TRAN))
    OK( GrB_Col_extract(x, GrB_NULL, GrB_NULL, M, GrB_ALL, nv(g), v, inp0_tran_desc) )
    nbrs, _ = findnz(x)
    OK( GrB_free(x) )
    OK( GrB_free(inp0_tran_desc) )
    return nbrs
end

function edges(g::BLASGraph)
    M = UpperTriangular(g.A)
    I, J, _ = findnz(M)
    OK( GrB_free(M) )
    return (SimpleEdge(I[i], J[i]) for i in 1:length(I))
end

function edges(g::BLASDiGraph)
    I, J, _ = findnz(g.A)
    return (SimpleEdge(I[i], J[i]) for i in 1:length(I))
end

zero(g::T) where T <: AbstractBLASGraph = T(GrB_Matrix([], [], []))

eltype(::AbstractBLASGraph) = Int64

edgetype(::AbstractBLASGraph) = SimpleEdge{Int64}

weights(g::AbstractBLASGraph) = BLASGraphWeights(g.A)

function dropzeros!(M::GrB_Matrix)
    outp_replace_desc = GrB_Descriptor(Dict(GrB_OUTP => GrB_REPLACE))
    OK( GrB_assign(M, M, GrB_NULL, M, GrB_ALL, 0, GrB_ALL, 0, outp_replace_desc) )
    OK( GrB_free(outp_replace_desc) )
end

# Algorithms
export BLASGraph, BLASDiGraph, free
export count_triangles
import LightGraphs.gdistances
include("algorithms/bfs.jl")
include("algorithms/tricount.jl")

end # module
