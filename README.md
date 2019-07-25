[![Build Status](https://travis-ci.org/abhinavmehndiratta/LightGraphsGraphBLAS.jl.svg?branch=master)](https://travis-ci.org/abhinavmehndiratta/LightGraphsGraphBLAS.jl)
## LightGraphsGraphBLAS.jl

### Examples:

```julia
julia> using GraphBLASInterface, SuiteSparseGraphBLAS, LightGraphsGraphBLAS, LightGraphs

julia> GrB_init(GrB_NONBLOCKING)
GrB_SUCCESS::GrB_Info = 0

julia> g = BLASGraph{Int64}(5)    # edge weights are of type Int64, the eltype for all graphs is UInt64 and cannot be changed
{5, 0} undirected graph

julia> add_edge!(g, 2, 4, 5)
true

julia> add_edge!(g, 1, 2, 7)
true

julia> foreach(println, edges(g))
Edge 1 => 2 with weight 7
Edge 2 => 4 with weight 5

julia> rem_edge!(g, 2, 4)
true

julia> foreach(println, edges(g))
Edge 1 => 2 with weight 7

julia> add_edge!(g, 1, 2, 9)    # reweight the edge
true

julia> foreach(println, edges(g))
Edge 1 => 2 with weight 9
```

**Create graph from AbstractSimpleGraph :**
```julia
julia> g = BLASDiGraph(SimpleDiGraph(5, 6))
{5, 6} directed graph

julia> foreach(println, edges(g))
Edge 1 => 3 with weight 1
Edge 1 => 4 with weight 1
Edge 1 => 5 with weight 1
Edge 2 => 5 with weight 1
Edge 3 => 5 with weight 1
Edge 4 => 2 with weight 1
```

**Create graph from GraphBLAS matrix :**
```julia
julia> I = OneBasedIndex[1, 2, 3]; J = OneBasedIndex[2, 3, 1]; X = Float64[7.2, 3.4, 5.6];

julia> M = GrB_Matrix(I, J, X)
GrB_Matrix{Float64}

julia> g = BLASDiGraph(M)
{3, 3} directed graph

julia> foreach(println, edges(g))
Edge 1 => 2 with weight 7.2
Edge 2 => 3 with weight 3.4
Edge 3 => 1 with weight 5.6

julia> outneighbors(g, 1)    # vertices of the graph are always of type UInt64 since they are indices of a GraphBLAS matrix
1-element Array{UInt64,1}:
 0x0000000000000002

julia> inneighbors(g, 1)
1-element Array{UInt64,1}:
 0x0000000000000003
```

**Create a graph from edge list :**
```julia
julia> e1 = SimpleWeightedEdge(1, 2, Int32(4))
Edge 1 => 2 with weight 4

julia> e2 = SimpleWeightedEdge(2, 4, Int32(8))
Edge 2 => 4 with weight 8

julia> edge_list = [e1, e2]
2-element Array{SimpleWeightedEdge{Int64,Int32},1}:
 Edge 1 => 2 with weight 4
 Edge 2 => 4 with weight 8

julia> g = BLASGraph(edge_list)
{4, 2} undirected graph

julia> get_weight(g, 1, 2)
4
```

**Create a graph from adjacency matrix :**
```julia
julia> A = [0 1 1; 1 0 0; 0 1 0]
3×3 Array{Int64,2}:
 0  1  1
 1  0  0
 0  1  0

julia> g = BLASDiGraph(A)
{3, 4} directed graph

julia> indegree(g, 1)
1

julia> outdegree(g, 1)
2

julia> using SparseArrays

julia> B = sparse(A)
3×3 SparseMatrixCSC{Int64,Int64} with 4 stored entries:
  [2, 1]  =  1
  [1, 2]  =  1
  [3, 2]  =  1
  [1, 3]  =  1

julia> g = BLASDiGraph(B)
{3, 4} directed graph
```
