function bellman_ford_shortest_paths(g::BLASGraph{T}, s::Union{Int64, UInt64}) where T
    s = OneBasedIndex(s)
    A = g.A
 
    # this implementation of the algorithm requires A[i, i] = 0 for all 0 <= i < n
    for i = 1:nv(g)
        u = OneBasedIndex(i)
        A[u, u] = zero(T)
    end

    # Initialize distance vector
    d = GrB_Vector(Int64, nv(g))
    d[s] = 0            # change d[s] to 0

    # tmp vector to store distance vector after n (i.e., V) loops
    dtmp = copy(d)

    iter = 0            # number of iterations
    same = false        # variable indicating if d = dtmp

    # terminate when no new path is found or more than n-1 loops
    while !same && iter < nv(g)-1
        # excute semiring on d and A, and save the result to d
        OK( GrB_vxm(dtmp, GrB_NULL, GrB_NULL, GxB_MIN_PLUS_INT64, d, A, GrB_NULL) )

        same = dtmp == d

        if !same
            dtmp, d = d, dtmp
        end

        iter += 1
    end

    # check for negative-weight cycle
    if !same
        OK( GrB_vxm(dtmp, GrB_NULL, GrB_NULL, GxB_MIN_PLUS_INT64, d, A, GrB_NULL) )
        if dtmp != d
            error("Negative-weight cycle found")
        end
    end

    dropzeros!(A)
    OK( GrB_free(dtmp) )
    return d
end
