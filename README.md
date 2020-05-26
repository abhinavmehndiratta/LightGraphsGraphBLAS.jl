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

On 6 cores -

**gdistances (algorithm: [bfs_simple](https://github.com/GraphBLAS/LAGraph/blob/master/Source/Algorithm/LAGraph_bfs_simple.c)) :**
```julia
julia> using GraphBLASInterface, SuiteSparseGraphBLAS, LightGraphsGraphBLAS, LightGraphs, MatrixDepot, BenchmarkTools, SNAPDatasets

julia> md = mdopen("DIMACS10/caidaRouterLevel");

julia> A = md.A

julia> lg = SimpleGraph(A)
{192244, 609066} undirected simple Int64 graph

julia> GrB_init(GrB_NONBLOCKING)
GrB_SUCCESS::GrB_Info = 0

julia> bg = BLASGraph(lg)
{192244, 609066} undirected graph

julia> @benchmark gdistances(lg, source) setup = (source = rand(1:nv(lg)))
BenchmarkTools.Trial: 
  memory estimate:  4.42 MiB
  allocs estimate:  9
  --------------
  minimum time:     159.524 μs (0.00% GC)
  median time:      29.321 ms (0.00% GC)
  mean time:        29.411 ms (0.00% GC)
  maximum time:     34.208 ms (0.00% GC)
  --------------
  samples:          170
  evals/sample:     1

julia> @benchmark gdistances(bg, source) setup = (source = rand(1:nv(bg)))
BenchmarkTools.Trial: 
  memory estimate:  2.55 KiB
  allocs estimate:  95
  --------------
  minimum time:     42.370 ms (0.00% GC)
  median time:      45.216 ms (0.00% GC)
  mean time:        45.386 ms (0.00% GC)
  maximum time:     53.582 ms (0.00% GC)
  --------------
  samples:          111
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
  minimum time:     9.511 ms (0.00% GC)
  median time:      10.535 ms (0.00% GC)
  mean time:        10.542 ms (0.20% GC)
  maximum time:     12.301 ms (0.00% GC)
  --------------
  samples:          473
  evals/sample:     1

julia> @benchmark gdistances(bg, source) setup = (source = rand(1:nv(bg)))
BenchmarkTools.Trial: 
  memory estimate:  1.56 KiB
  allocs estimate:  60
  --------------
  minimum time:     16.929 ms (0.00% GC)
  median time:      20.320 ms (0.00% GC)
  mean time:        20.280 ms (0.00% GC)
  maximum time:     35.234 ms (0.00% GC)
  --------------
  samples:          247
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
  minimum time:     425.521 μs (0.00% GC)
  median time:      429.602 μs (0.00% GC)
  mean time:        442.115 μs (0.29% GC)
  maximum time:     1.709 ms (62.72% GC)
  --------------
  samples:          10000
  evals/sample:     1

julia> @benchmark gdistances(bg, source) setup = (source = rand(1:nv(bg)))
BenchmarkTools.Trial: 
  memory estimate:  1.13 KiB
  allocs estimate:  44
  --------------
  minimum time:     662.463 μs (0.00% GC)
  median time:      1.472 ms (0.00% GC)
  mean time:        1.478 ms (0.00% GC)
  maximum time:     10.979 ms (0.00% GC)
  --------------
  samples:          3313
  evals/sample:     1
```
