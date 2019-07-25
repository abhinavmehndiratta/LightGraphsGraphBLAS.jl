# check if GraphBLAS operation was successful, else throw an error
function OK(info::GrB_Info)
    info != GrB_SUCCESS && error(info)
    return true
end
