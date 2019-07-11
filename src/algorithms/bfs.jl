function gdistances(g::BLASGraph, s::Int64)
    A = g.A
    desc = GrB_Descriptor(Dict(GrB_MASK => GrB_SCMP, GrB_OUTP => GrB_REPLACE))
    n = nv(g)

    v = GrB_Vector(Int64, n)        # result vector
    OK( GrB_assign(v, GrB_NULL, GrB_NULL, 0, GrB_ALL, n, GrB_NULL) )

    q = GrB_Vector(Bool, n)         # nodes visited at each level
    q[s] = true

    successor = true        # true when some successor found
    level = 1
    while successor && level <= n
        # v<q> = level, using vector assign with q as the mask
        OK( GrB_assign(v, q, GrB_NULL, level, GrB_ALL, n, GrB_NULL) )

        # q<!v> = q ||.&& A ; finds all the unvisited
        # successors from current q, using !v as the mask
        OK( GrB_vxm(q, v, GrB_NULL, GxB_LOR_LAND_BOOL, q, A, desc) )

        # successor = ||(q)
        successor = GrB_reduce(GxB_LOR_BOOL_MONOID, q, GrB_NULL)

        level += 1
    end

    subtract1(x) = x-1
    sub1 = GrB_UnaryOp()
    GrB_UnaryOp_new(sub1, subtract1, GrB_INT64, GrB_INT64)

    OK( GrB_apply(v, v, GrB_NULL, sub1, v, GrB_NULL) )

    mask_scmp_desc = GrB_Descriptor(Dict(GrB_MASK => GrB_SCMP))
    OK( GrB_assign(v, v, GrB_NULL, typemax(Int64), GrB_ALL, 0, mask_scmp_desc) )
    
    OK( GrB_free(q) )
    OK( GrB_free(desc) )
    OK( GrB_free(sub1) )
    OK( GrB_free(mask_scmp_desc) )

    v[s] = 0

    return v
end
