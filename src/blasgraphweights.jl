mutable struct BLASGraphWeights{T} <: AbstractMatrix{T}
    A::GrB_Matrix{T}
end
Base.show(::IO, ::MIME{Symbol("text/plain")}, ::BLASGraphWeights) = print("GraphBLAS weight matrix")
size(w::BLASGraphWeights) = size(w.A)
getindex(w::BLASGraphWeights, r::Integer, c::Integer) = getindex(w.A, OneBasedIndex(r), OneBasedIndex(c))
