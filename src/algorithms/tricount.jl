function count_triangles(g::BLASGraph)
    M = g.A
    L = LowerTriangular(M)
    C = GrB_Matrix(Int64, size(M)...)

    # Descriptor for mxm
    desc = GrB_Descriptor(GrB_INP1 => GrB_TRAN) # transpose the second matrix

    OK( GrB_mxm(C, L, GrB_NULL, GxB_PLUS_TIMES_INT64, L, L, desc) ) # C<L> = L ∗.+ L'
    ntriangles = GrB_reduce(GxB_PLUS_INT64_MONOID, C, GrB_NULL)

    GrB_free(C)
    GrB_free(L)
    GrB_free(desc)

    return ntriangles
end
