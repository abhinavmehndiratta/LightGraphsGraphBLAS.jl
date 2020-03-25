using LightGraphs, GraphBLASInterface, SuiteSparseGraphBLAS, LightGraphsGraphBLAS, SparseArrays, Test

GrB_init(GrB_NONBLOCKING)
I = OneBasedIndex[1, 2, 3]; J = OneBasedIndex[2, 3, 1]; X = Float64[7.2, 3.4, 5.6]
M = GrB_Matrix(I, J, X)
g = BLASDiGraph(M)
e1 = SimpleWeightedEdge(1, 2, 7.2)
e2 = SimpleWeightedEdge(2, 3, 3.4)
e3 = SimpleWeightedEdge(3, 1, 5.6)
edge_list = [e1, e2, e3]
@test edges(g) == edge_list
