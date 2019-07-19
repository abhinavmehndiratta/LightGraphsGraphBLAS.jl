import Base:
    size, getindex

struct BLASGraphWeights{T} <: AbstractMatrix{T}
    A::GrB_Matrix{T}
end

size(w::BLASGraphWeights) = size(w.A)
getindex(w::BLASGraphWeights, r::Union{Int64, UInt64}, c::Union{Int64, UInt64}) = getindex(w.A, OneBasedIndex(r), OneBasedIndex(c))
