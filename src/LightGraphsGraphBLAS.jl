module LightGraphsGraphBLAS

using LightGraphs, GraphBLASInterface, SuiteSparseGraphBLAS, SparseArrays, LinearAlgebra

import LightGraphs:
    ne, nv, is_directed, has_vertex, vertices, has_edge, add_edge!, rem_edge!,
    indegree, outdegree, inneighbors, outneighbors, edgetype, edges, weights,
    gdistances, bellman_ford_shortest_paths,prim_mst

import SimpleWeightedGraphs:
    SimpleWeightedEdge, weight

import Base:
        show, size, getindex, eltype, copy, argmin

include("abstractblasgraph.jl")
include("blasgraph.jl")
include("blasdigraph.jl")
include("blasgraphweights.jl")
include("algorithms/bfs.jl")
include("algorithms/tricount.jl")
include("algorithms/bellman_ford.jl")
include("algorithms/prim_mst.jl")
include("utils.jl")

export BLASGraph, BLASDiGraph, SimpleWeightedEdge, get_weight

export count_triangles

end # module
