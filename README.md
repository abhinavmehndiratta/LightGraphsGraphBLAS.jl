[![Build Status](https://travis-ci.org/abhinavmehndiratta/LightGraphsGraphBLAS.jl.svg?branch=master)](https://travis-ci.org/abhinavmehndiratta/LightGraphsGraphBLAS.jl)
## LightGraphsGraphBLAS.jl

The edge weights can be of type `Bool, Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32 or Float64` (i.e., the GraphBLAS predefined types). User-defined types are not supported.

### Examples

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

### Benchmarks

**gdistances (algorithm: [bfs_simple](https://github.com/GraphBLAS/LAGraph/blob/master/Source/Algorithm/LAGraph_bfs_simple.c)) :**
```julia
julia> using GraphBLASInterface, SuiteSparseGraphBLAS, LightGraphsGraphBLAS, LightGraphs, MatrixDepot, BenchmarkTools, SNAPDatasets

julia> md = mdopen("DIMACS10/caidaRouterLevel");

julia> A = md.A

julia> lg = SimpleGraph(A)
{192244, 609066} undirected simple Int64 graph

julia> bg = BLASGraph(lg)
{192244, 609066} undirected graph

julia> @benchmark gdistances(lg, source) setup = (source = rand(1:nv(lg)))
BenchmarkTools.Trial:
  memory estimate:  4.42 MiB
  allocs estimate:  11
  --------------
  minimum time:     212.493 μs (0.00% GC)
  median time:      39.958 ms (0.00% GC)
  mean time:        40.258 ms (0.00% GC)
  maximum time:     45.352 ms (0.00% GC)
  --------------
  samples:          124
  evals/sample:     1

julia> @benchmark gdistances(bg, source) setup = (source = rand(1:nv(bg)))
BenchmarkTools.Trial:
  memory estimate:  2.91 KiB
  allocs estimate:  96
  --------------
  minimum time:     85.847 ms (0.00% GC)
  median time:      87.728 ms (0.00% GC)
  mean time:        88.467 ms (0.00% GC)
  maximum time:     94.347 ms (0.00% GC)
  --------------
  samples:          57
  evals/sample:     1

julia> lg = loadsnap(:soc_slashdot0902_u)
{82168, 582533} undirected simple Int64 graph

julia> bg = BLASGraph(lg)
{82168, 582533} undirected graph

julia> @benchmark gdistances(lg, source) setup = (source = rand(1:nv(lg)))
BenchmarkTools.Trial:
  memory estimate:  1.89 MiB
  allocs estimate:  8
  --------------
  minimum time:     13.377 ms (0.00% GC)
  median time:      15.876 ms (0.00% GC)
  mean time:        15.944 ms (0.25% GC)
  maximum time:     21.372 ms (0.00% GC)
  --------------
  samples:          313
  evals/sample:     1

julia> @benchmark gdistances(bg, source) setup = (source = rand(1:nv(bg)))
BenchmarkTools.Trial:
  memory estimate:  2.06 KiB
  allocs estimate:  66
  --------------
  minimum time:     36.501 ms (0.00% GC)
  median time:      39.880 ms (0.00% GC)
  mean time:        40.116 ms (0.00% GC)
  maximum time:     53.307 ms (0.00% GC)
  --------------
  samples:          125
  evals/sample:     1

julia> lg = loadsnap(:facebook_combined)
{4039, 88234} undirected simple Int64 graph

julia> bg = BLASGraph(lg)
{4039, 88234} undirected graph

julia> @benchmark gdistances(lg, source) setup = (source = rand(1:nv(lg)))
BenchmarkTools.Trial:
  memory estimate:  95.69 KiB
  allocs estimate:  8
  --------------
  minimum time:     634.694 μs (0.00% GC)
  median time:      683.736 μs (0.00% GC)
  mean time:        693.038 μs (0.30% GC)
  maximum time:     1.699 ms (53.43% GC)
  --------------
  samples:          6953
  evals/sample:     1

julia> @benchmark gdistances(bg, source) setup = (source = rand(1:nv(bg)))
BenchmarkTools.Trial:
  memory estimate:  1.63 KiB
  allocs estimate:  50
  --------------
  minimum time:     2.114 ms (0.00% GC)
  median time:      2.462 ms (0.00% GC)
  mean time:        2.462 ms (0.00% GC)
  maximum time:     4.275 ms (0.00% GC)
  --------------
  samples:          2018
  evals/sample:     1
```
