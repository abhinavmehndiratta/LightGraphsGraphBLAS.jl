# check if GraphBLAS operation was successful, else throw an error
function OK(info::GrB_Info)
    info != GrB_SUCCESS && error(info)
    return true
end

function PlusOp(T::DataType)
    if T == Bool
        return GrB_PLUS_BOOL
    elseif T == Int8
        return GrB_PLUS_INT8
    elseif T == UInt8
        return GrB_PLUS_UINT8
    elseif T == Int16
        return GrB_PLUS_INT16
    elseif T == UInt16
        return GrB_PLUS_UINT16
    elseif T == Int32
        return GrB_PLUS_INT32
    elseif T == UInt32
        return GrB_PLUS_UINT32
    elseif T == Int64
        return GrB_PLUS_INT64
    elseif T == UInt64
        return GrB_PLUS_UINT64
    elseif  T == Float32
        return GrB_PLUS_FP32
    end
    return GrB_PLUS_FP64
end

function MinOp(T::DataType)
    if T == Bool
        return GrB_MIN_BOOL
    elseif T == Int8
        return GrB_MIN_INT8
    elseif T == UInt8
        return GrB_MIN_UINT8
    elseif T == Int16
        return GrB_MIN_INT16
    elseif T == UInt16
        return GrB_MIN_UINT16
    elseif T == Int32
        return GrB_MIN_INT32
    elseif T == UInt32
        return GrB_MIN_UINT32
    elseif T == Int64
        return GrB_MIN_INT64
    elseif T == UInt64
        return GrB_MIN_UINT64
    elseif  T == Float32
        return GrB_MIN_FP32
    end
    return GrB_MIN_FP64
end

function IdentityOp(T::DataType)
    if T == Bool
        return GrB_IDENTITY_BOOL
    elseif T == Int8
        return GrB_IDENTITY_INT8
    elseif T == UInt8
        return GrB_IDENTITY_UINT8
    elseif T == Int16
        return GrB_IDENTITY_INT16
    elseif T == UInt16
        return GrB_IDENTITY_UINT16
    elseif T == Int32
        return GrB_IDENTITY_INT32
    elseif T == UInt32
        return GrB_IDENTITY_UINT32
    elseif T == Int64
        return GrB_IDENTITY_INT64
    elseif T == UInt64
        return GrB_IDENTITY_UINT64
    elseif  T == Float32
        return GrB_IDENTITY_FP32
    end
    return GrB_IDENTITY_FP64
end
