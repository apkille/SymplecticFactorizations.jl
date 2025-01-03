struct Takagi{T,M<:AbstractArray{T},N<:AbstractVector{<:Real}} <: Factorization{T}
    Q::M
    S::N
    function Takagi{T,M,N}(Q, S) where {T,M<:AbstractArray{T},N<:AbstractVector{<:Real}}
        require_one_based_indexing(Q, S)
        new{T,M,N}(Q, S)
    end
end
Takagi{T}(Q::AbstractArray{T}, S::AbstractVector{<:Real}) where {T} = Takagi{T,typeof(Q),typeof(S)}(Q,S)

# iteration for destructuring into components
Base.iterate(F::Takagi) = (F.Q, Val(:S))
Base.iterate(F::Takagi, ::Val{:S}) = (F.S, Val(:done))
Base.iterate(F::Takagi, ::Val{:done}) = nothing

function takagi(x::AbstractMatrix{T}) where {T<:Union{Real,Complex}}
    fact = svd(Symmetric(x))
    Tt = T <: Complex ? T : ComplexF64
    Q = zeros(Tt, size(x))
    mul!(Q, fact.U, sqrt(conj(transpose(fact.U) * fact.V)))
    return Takagi{Tt}(Q, fact.S)
end
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::Takagi{<:Any,<:AbstractArray,<:AbstractVector})
    Base.summary(io, F); println(io)
    println(io, "Q factor:")
    Base.show(io, mime, F.Q)
    println(io, "\nsingular values:")
    Base.show(io, mime, F.S)
end