struct BLASGraphWeights{T} <: AbstractMatrix{T}
    A::GrB_Matrix{T}
end
show(io::IO, ::BLASGraphWeights) = print("Weight matrix")
size(w::BLASGraphWeights) = size(w.A)
getindex(w::BLASGraphWeights, r::Union{Int64, UInt64}, c::Union{Int64, UInt64}) = getindex(w.A, OneBasedIndex(r), OneBasedIndex(c))
