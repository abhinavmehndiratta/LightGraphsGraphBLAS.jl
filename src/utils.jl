# check if GraphBLAS operation was successful, else throw an error
function OK(info::GrB_Info)
    info != GrB_SUCCESS && error(info)
    return true
end

struct Dists
    v::GrB_Vector
end

function show(io::IO, w::Dists)
    n = Int64(size(w.v, 1))
    println(n, "-element GraphBLAS distance vector:")
    if n <= 20
        for i = 1:n-1
            println(i, " => ", w[i])
        end
        print(n, " => ", w[n])
    else
        for i = 1:10
            println(i, " => ", w[i])
        end
        println("⋮")
        for i = (n-10+1):n-1
            println(i, " => ", w[i])
        end
        print(n, " => ", w[n])
    end    
end

getindex(w::Dists, i::Union{Int64, UInt64}) = getindex(w.v, OneBasedIndex(i))
