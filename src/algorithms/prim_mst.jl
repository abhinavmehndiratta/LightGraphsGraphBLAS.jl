# TO-DO can be optimized when user-defined types are supported
function argmin(v::GrB_Vector)
	n = SuiteSparseGraphBLAS.nnz(v)
    n == 0 && error(GrB_PANIC)
    idx, vals = SuiteSparseGraphBLAS.findnz(v)
    return idx[argmin(vals)]
end

# return sum of edge-weights of MST using Prim's algorithm
function prim_mst(graph::BLASGraph{T}) where T
	A = graph.A
    START_NODE = ZeroBasedIndex(0)
    weight = zero(T)
	nvg = nv(graph)
	desc1 = GrB_Descriptor(Dict(GrB_MASK => GrB_SCMP, GrB_OUTP => GrB_REPLACE))

    mask = GrB_Vector(Bool, nvg)
    mask[START_NODE] = true

    s = GrB_Vector(T, nvg)
    OK( GrB_assign(s, mask, GrB_NULL, 0, GrB_ALL, 0, desc1) )

    d = GrB_Vector(T, nvg)
    OK( GrB_extract(d, GrB_NULL, GrB_NULL, A, GrB_ALL, 0, START_NODE, GrB_NULL) )

    temp = GrB_Vector(T, nvg)
    Arow = GrB_Vector(T, nvg)

    PlusT = PlusOp(T); IdentityT = IdentityOp(T); MinT = MinOp(T)

    while SuiteSparseGraphBLAS.nnz(mask) < nvg
        empty!(temp)

        OK( GrB_eWiseMult(temp, GrB_NULL, GrB_NULL, PlusT, s, d, GrB_NULL) )

        u = argmin(temp)
        weight += d[u]
        mask[u] = true

        OK( GrB_apply(s, mask, GrB_NULL, IdentityT, s, desc1) )

        empty!(Arow)

        OK( GrB_extract(Arow, GrB_NULL, GrB_NULL, A, GrB_ALL, 0, u, GrB_NULL) )

        OK( GrB_eWiseAdd(d, GrB_NULL, GrB_NULL, MinT, d, Arow, GrB_NULL) )
    end

    OK( GrB_free(desc1) )
    OK( GrB_free(mask) )
    OK( GrB_free(s) )
    OK( GrB_free(d) )
    OK( GrB_free(temp) )
    OK( GrB_free(Arow) )
    return weight
end
