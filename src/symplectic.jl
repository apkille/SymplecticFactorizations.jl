struct Symplectic{F<:SymplecticForm, T<:Number, D<:AbstractMatrix{<:T}} <: AbstractMatrix{T}
    form::F
    data::D
    function Symplectic{F,T,D}(form, data) where {F<:SymplecticForm,T<:Number,D<:AbstractMatrix{<:T}}
        LinearAlgebra.require_one_based_indexing(data)
        new{F,T,D}(form, data)
    end
end

"""
    Symplectic(J::SymplecticForm, S::AbstractMatrix) <: AbstractMatrix
    Symplectic(J::SymplecticForm, S::AbstractMatrix; atol = 1e-8, rtol = atol) <: AbstractMatrix

Construct a wrapper of a symplectic matrix `S` with its corresponding symplectic basis
defined by the symplectic form `J`.
"""
Symplectic(form::F, data::D) where {F<:SymplecticForm, D<:AbstractMatrix} = Symplectic{F,eltype(D),D}(form, data)
Symplectic(x::Symplectic) = x

Base.size(x::Symplectic) = size(x.data)
Base.size(x::Symplectic, n) = size(x.data, n)
Base.axes(x::Symplectic) = axes(x.data)
Base.eltype(x::Symplectic) = eltype(x.data)

Base.@propagate_inbounds Base.getindex(x::Symplectic, i::Int) = x.data[i]
Base.@propagate_inbounds Base.getindex(x::Symplectic, i::Int, j::Int) = x.data[i,j]
Base.@propagate_inbounds Base.getindex(x::Symplectic, I::Vararg{Int,N}) where {N} = getindex(x.data, I)
Base.@propagate_inbounds Base.setindex!(x::Symplectic, v, i::Int) = setindex!(x.data, v, i)
Base.@propagate_inbounds Base.setindex!(x::Symplectic, v, I::Vararg{Int, N}) where {N} = setindex!(x, v, I)

Base.similar(x::Symplectic, ::Type{T}) where {T} = similar(x.data, T)
Base.similar(x::Symplectic, dims::Dims{N}) where {N} = similar(x.data, dims)
Base.similar(x::Symplectic, ::Type{T}, dims::Dims{N}) where {T,N} = similar(x.data, T, dims)

Base.Matrix(x::Symplectic) = Matrix(x.data)
Base.Array(x::Symplectic) = Matrix(x)
Base.AbstractMatrix{T}(x::Symplectic) where {T} = Symplectic(x.form, AbstractMatrix{T}(x.data))
Base.parent(x::Symplectic) = x.data

function issymplectic(form::SymplecticForm, x::T; atol::R1 = 0, rtol::R2 = atol) where {T,R1<:Real,R2<:Real}
    omega = symplecticform(form)
    return isapprox(x * omega * x', omega; atol = atol, rtol = rtol)
end

"""
    randsymplectic(form::SymplecticForm)

Calculate a random symplectic matrix in symplectic representation defined by `basis`.
"""
randsymplectic(form::SymplecticForm{N}) where {N<:Int} = _randsymplectic(form)
function _randsymplectic(form::PairForm{N}) where {N<:Int}
    n = form.n
    # random Block-Messiah/Euler decomposition
    O, O′ = _rand_orthogonal_symplectic(form), _rand_orthogonal_symplectic(form)
    rs = rand(n)
    D = Diagonal(collect(Iterators.flatten((i, 1/i) for i in rs)))
    return O * D * O′
end
function _randsymplectic(form::BlockForm{N}) where {N<:Int}
    n = form.n
    # random Block-Messiah/Euler decomposition
    O, O′ = _rand_orthogonal_symplectic(form), _rand_orthogonal_symplectic(form)
    rs = rand(n)
    D = Diagonal(vcat(rs, 1 ./ rs))
    return O * D * O′
end

# Generates random orthogonal symplectic matrix by blocking real
# and imaginary parts of a random unitary matrix
function _rand_orthogonal_symplectic(form::PairForm{N}) where {N<:Int}
    n = form.n
    U = _rand_unitary(form)
    O = zeros(2*n, 2*n)
    @inbounds for i in Base.OneTo(n), j in Base.OneTo(n)
        val = U[i,j]
        O[2*i-1, 2*j-1] = real(val)
        O[2*i, 2*j-1] = -imag(val)
        O[2*i-1, 2*j] = imag(val)
        O[2*i, 2*j] = real(val)
    end
    return O
end
function _rand_orthogonal_symplectic(form::BlockForm{N}) where {N<:Int}
    n = form.n
    U = _rand_unitary(form)
    O = zeros(2*n, 2*n)
    @inbounds for i in Base.OneTo(n), j in Base.OneTo(n)
        val = U[i,j]
        O[i, j] = real(val)
        O[i+n, j] = -imag(val)
        O[i, j+n] = imag(val)
        O[i+n, j+n] = real(val)
    end
    return O
end

# Generates unitary matrix randomly distributed over Haar measure;
# see https://arxiv.org/abs/math-ph/0609050 for algorithm.
# This approach is faster and creates less allocations than rand(Haar(2), n) from RandomMatrices.jl
function _rand_unitary(form::SymplecticForm{N}) where {N<:Int}
    n = form.n
    M = rand(ComplexF64, n, n) ./ sqrt(2.0)
    q, r = qr(M)
    d = Diagonal([r[i, i] / abs(r[i, i]) for i in Base.OneTo(n)])
    return q * d
end