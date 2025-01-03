struct Williamson{T,M<:AbstractArray{T},N<:AbstractVector{T}} <: Factorization{T}
    S::M
    spectrum::N
    function Williamson{T,M,N}(S, spectrum) where {T,M<:AbstractArray{T},N<:AbstractVector{T}}
        require_one_based_indexing(S, spectrum)
        new{T,M,N}(S, spectrum)
    end
end
Williamson{T}(S::AbstractArray{T}, spectrum::AbstractVector{T},) where {T} = Williamson{T,typeof(S),typeof(spectrum)}(S, spectrum)

# iteration for destructuring into components
Base.iterate(F::Williamson) = (F.S, Val(:spectrum))
Base.iterate(F::Williamson, ::Val{:spectrum}) = (F.spectrum, Val(:done))
Base.iterate(F::Williamson, ::Val{:done}) = nothing

function williamson(form::PairForm, x::AbstractMatrix{T}) where {T<:Real}
    J = symplecticform(form)
    spectrum = filter(i -> i > 0, imag.(eigvals(J * x, sortby = λ -> abs(λ))))
    D = Diagonal(repeat(spectrum, inner = 2))
    sqrtV = Symmetric(x)^(-1//2)
    X = sqrtV * J * sqrtV
    U = eigvecs(X, sortby = λ -> -abs(λ))
    G = zeros(ComplexF64, 2*form.n, 2*form.n)
    @inbounds for i in Base.OneTo(form.n)
        G[2i-1, 2i-1] = -im/sqrt(2)
        G[2i-1, 2i] = im/sqrt(2)
        G[2i, 2i-1] = 1/sqrt(2)
        G[2i, 2i] = 1/sqrt(2)
    end
    R = real(G * U')
    S = sqrt(D) * R * sqrtV
    return Williamson{T}(S, spectrum)
end
function williamson(form::BlockForm, x::AbstractMatrix{T}) where {T<:Real}
    J = symplecticform(form)
    spectrum = filter(i -> i > 0, imag.(eigvals(J * x, sortby = λ -> abs(λ))))
    D = Diagonal(repeat(spectrum, 2))
    sqrtV = Symmetric(x)^(-1//2)
    X = sqrtV * J * sqrtV
    U = eigvecs(X, sortby = λ -> (-sign(imag(λ)), -abs(λ)))
    G = zeros(ComplexF64, 2*form.n, 2*form.n)
    @inbounds for i in Base.OneTo(form.n)
        G[i, i] = -im/sqrt(2)
        G[i, i+form.n] = im/sqrt(2)
        G[i+form.n, i] = 1/sqrt(2)
        G[i+form.n, i+form.n] = 1/sqrt(2)
    end
    R = real(G * U')
    S = sqrt(D) * R * sqrtV
    return Williamson{T}(S, spectrum)
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, F::Williamson{<:Any,<:AbstractArray,<:AbstractVector})
    Base.summary(io, F); println(io)
    println(io, "S factor:")
    Base.show(io, mime, F.S)
    println(io, "\nsymplectic spectrum:")
    Base.show(io, mime, F.spectrum)
end