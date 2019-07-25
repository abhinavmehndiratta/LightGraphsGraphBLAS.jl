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
            println(w[i])
        end
        print(w[n])
    else
        for i = 1:10
            println(w[i])
        end
        println("⋮")
        for i = (n-10+1):n-1
            println(w[i])
        end
        print(w[n])
    end    
end

getindex(w::Dists, i::Integer) = getindex(w.v, OneBasedIndex(i))