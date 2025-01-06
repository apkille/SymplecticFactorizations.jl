struct Polar{T,M<:AbstractArray{T},N<:AbstractArray{T}} <: Factorization{T}
    O::M
    P::N
    function Polar{T,M,N}(O, P) where {T,M<:AbstractArray{T},N<:AbstractArray{T}}
        require_one_based_indexing(O, P)
        new{T,M,N}(O, P)
    end
end
Polar{T}(O::AbstractArray{T}, P::AbstractArray{T}) where {T} = Polar{T,typeof(O),typeof(P)}(O,P)

# iteration for destructuring into components
Base.iterate(F::Polar) = (F.O, Val(:P))
Base.iterate(F::Polar, ::Val{:P}) = (F.P, Val(:done))
Base.iterate(F::Polar, ::Val{:done}) = nothing

function polar(x::AbstractMatrix{T}) where {T}
    fact = svd(x)
    dims = size(x)
    O, P = zeros(T, dims), zeros(T, dims)
    mul!(O, fact.U, fact.Vt)
    copyto!(P, fact.V * Diagonal(fact.S) * fact.Vt)
    return Polar{T}(O, P)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::Polar{<:Any,<:AbstractArray,<:AbstractArray})
    Base.summary(io, F); println(io)
    println(io, "O factor:")
    Base.show(io, mime, F.O)
    println(io, "\nP factor:")
    Base.show(io, mime, F.P)
end