module LightGraphsGraphBLAS

using LightGraphs, GraphBLASInterface, SuiteSparseGraphBLAS

import LightGraphs:
    ne, nv, is_directed, has_vertex, vertices, has_edge, add_edge!, rem_edge!,
    indegree, outdegree, inneighbors, outneighbors, edgetype, edges,
    SimpleGraphs.SimpleEdge, weights, gdistances, bellman_ford_shortest_paths

import Base.eltype

include("abstractblasgraph.jl")
include("blasgraph.jl")
include("blasdigraph.jl")
include("blasgraphweights.jl")
include("algorithms/bfs.jl")
include("algorithms/tricount.jl")
include("algorithms/bellman_ford.jl")
include("utils.jl")

export BLASGraph, BLASDiGraph

export count_triangles

end # module
