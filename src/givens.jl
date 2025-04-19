struct SymplecticGivens{F<:SymplecticForm,N<:Int,T} <: AbstractRotation{T}
    form::F
    k::N
    c::T
    s::T
    SymplecticGivens(form::F, k::N, c::T, s::T) where {F<:SymplecticForm,N<:Int,T} = new{F,N,T}(form, k, c, s)
end
SymplecticGivens(x::SymplecticGivens) = x

function givens(form::SymplecticForm, k::Int, θ::T) where {T}
    k < form.n || throw(ArgumentError("k must be less than n."))
    c, s = cos(θ), sin(θ)
    return SymplecticGivens(form, k, c, s)
end

function Base.Matrix(x::SymplecticGivens{F,N,T}) where {F<:BlockForm,N<:Int,T}
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    M = zeros(T, 2n, 2n)
    @inbounds for i in Base.OneTo(n)
        if i == k
            M[i,i] = c
            M[i+n,i+n] = c
            M[i,i+n] = s
            M[i+n,i] = -s
        else
            M[i,i] = oneunit(c)
            M[i+n,i+n] = oneunit(c)
        end
    end
    return M
end
function Base.Matrix(x::SymplecticGivens{F,N,T}) where {F<:PairForm,N<:Int,T}
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    M = Matrix{T}(I, 2n, 2n)
    M[2k-1,2k-1] = c
    M[2k-1,2k] = s
    M[2k,2k-1] = -s
    M[2k,2k] = c
    return M
end
Base.Array(x::SymplecticGivens) = Matrix(x)
function Base.AbstractMatrix{T1}(x::SymplecticGivens{F,N,T2}) where {F<:BlockForm,N<:Int,T1,T2}
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    M = zeros(T1, 2n, 2n)
    @inbounds for i in Base.OneTo(n)
        if i == k
            M[i,i] = c
            M[i+n,i+n] = c
            M[i,i+n] = s
            M[i+n,i] = -s
        else
            M[i,i] = oneunit(c)
            M[i+n,i+n] = oneunit(c)
        end
    end
    return M
end
function Base.AbstractMatrix{T1}(x::SymplecticGivens{F,N,T2}) where {F<:PairForm,N<:Int,T1,T2}
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    M = Matrix{T1}(I, 2n, 2n)
    M[2k-1,2k-1] = c
    M[2k-1,2k] = s
    M[2k,2k-1] = -s
    M[2k,2k] = c
    return M
end
Base.checkbounds(x::SymplecticGivens, i::Int, j::Int) = (form = x.form; i <= 2 * form.n && j <= 2 * form.n)
Base.checkbounds(x::SymplecticGivens, i::Int) = (form = x.form; i <= (2 * form.n)^2)
Base.isequal(x::SymplecticGivens, y::SymplecticGivens) = x.form == y.form && x.k == y.k && x.c == y.c && x.s == y.s
Base.isapprox(x::SymplecticGivens, y::SymplecticGivens) = x.form == y.form && isapprox(x.k, y.k) && isapprox(x.c, y.c) && isapprox(x.s, y.s)
Base.size(x::SymplecticGivens) = (form = x.form; (2 * form.n, 2 * form.n))
Base.size(x::SymplecticGivens, dim::Int) = (form = x.form; 2 * form.n)
Base.length(x::SymplecticGivens) = prod(size(x))
Base.axes(x::SymplecticGivens) = (form = x.form; (Base.OneTo(2 * form.n), Base.OneTo(2 * form.n)))
Base.eltype(x::SymplecticGivens) = eltype(x.c)

Base.@propagate_inbounds function Base.getindex(x::SymplecticGivens{F,N,T}, i::Int, j::Int) where {F<:BlockForm,N<:Int,T}
    @boundscheck checkbounds(x, i, j)
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    if i == j
        if i == k || i == (n+k)
            return c
        else
            return oneunit(c)
        end
    elseif i == k && j == (n+k)
        return s
    elseif i == (n+k) && j == k
        return -s
    else 
        return zero(c)
    end
end
Base.@propagate_inbounds function Base.getindex(x::SymplecticGivens{F,N,T}, i::Int) where {F<:BlockForm,N<:Int,T}
    @boundscheck checkbounds(x, i)
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    if i % (2n+1) == 1
        if (i-k) % (2n) == 0 || (i-n-k) % (2n) == 0
            return c
        else
            return oneunit(c)
        end
    elseif i == (2*n*k-n+k)
        return -s
    elseif i == (2*n^2+2n*(k-1)+k)
        return s
    else
        return zero(c)
    end
end
Base.@propagate_inbounds function Base.getindex(x::SymplecticGivens{F,N,T}, i::Int, j::Int) where {F<:PairForm,N<:Int,T}
    @boundscheck checkbounds(x, i, j)
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    if i == j
        if i == (2k-1) || i == (2k)
            return c
        else
            return oneunit(c)
        end
    elseif i == (2k-1) && j == (2k)
        return s
    elseif i == (2k) && j == (2k-1)
        return -s
    else 
        return zero(c)
    end
end
Base.@propagate_inbounds function Base.getindex(x::SymplecticGivens{F,N,T}, i::Int) where {F<:PairForm,N<:Int,T}
    @boundscheck checkbounds(x, i)
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    if i % (2n+1) == 1
        if i == (4*n*(k-1)+2k-1) || i == (4*n*(k-1)+2(n+k))
            return c
        else
            return oneunit(c)
        end
    elseif i == (4*n*(k-1)+2k)
        return -s
    elseif i == (4*n*(k-1)+2(n+k)-1)
        return s
    else
        return zero(c)
    end
end

LinearAlgebra.adjoint(x::SymplecticGivens) = SymplecticGivens(x.form, x.k, x.c, -x.s)
LinearAlgebra.inv(x::SymplecticGivens) = adjoint(x)
Base.copy(x::SymplecticGivens) = SymplecticGivens(x.form, copy(x.k), copy(x.c), copy(x.s))
@inline function Base.copyto!(dest::AbstractMatrix, src::SymplecticGivens{F,N,T}) where {F<:BlockForm,N<:Int,T}
    LinearAlgebra.require_one_based_indexing(dest)
    size(dest, 1) == size(dest, 2) || throw(ArgumentError("cannot copy a SymplecticGivens object to a non-square matrix."))
    n, k, c, s = (src.form).n, src.k, src.c, src.s
    @inbounds for i in Base.OneTo(n)
        if i == k
            dest[i,i] = c
            dest[i+n,i+n] = c
            dest[i,i+n] = s
            dest[i+n,i] = -s
        else
            dest[i,i] = oneunit(c)
            dest[i+n,i+n] = oneunit(c)
        end
    end
end
@inline function Base.copyto!(dest::AbstractMatrix, src::SymplecticGivens{F,N,T}) where {F<:PairForm,N<:Int,T}
    LinearAlgebra.require_one_based_indexing(dest)
    size(dest, 1) == size(dest, 2) || throw(ArgumentError("cannot copy a SymplecticGivens object to a non-square matrix."))
    n, k, c, s = (src.form).n, src.k, src.c, src.s
    @inbounds for i in Base.OneTo(n)
        if i == k
            dest[2i-1,2i-1] = c
            dest[2i,2i] = c
            dest[2i-1,2i] = s
            dest[2i,2i-1] = -s
        else
            dest[2i-1,2i-1] = oneunit(c)
            dest[2i,2i] = oneunit(c)
        end
    end
end

@inline function LinearAlgebra.lmul!(x::SymplecticGivens{F,N,T}, y::AbstractVecOrMat) where {F<:BlockForm,N<:Int,T}
    LinearAlgebra.require_one_based_indexing(y)
    size(y, 1) == size(y, 2) || throw(ArgumentError("cannot compute the matrix product between a SymplecticGivens object and a non-square matrix."))
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    @inbounds for i in Base.OneTo(n)
        y1, y2 = y[k,i], y[n+k,i]
        y1n, y2n = y[k,i+n], y[n+k,i+n]
        y[k,i] = c * y1 + s * y2
        y[k,i+n] = c * y1n + s * y2n
        y[n+k,i] = -s * y1 + c * y2
        y[n+k,i+n] = -s * y1n + c * y2n
    end
    return y
end
@inline function LinearAlgebra.rmul!(x::AbstractVecOrMat, y::SymplecticGivens{F,N,T}) where {F<:BlockForm,N<:Int,T}
    LinearAlgebra.require_one_based_indexing(x)
    size(x, 1) == size(x, 2) || throw(ArgumentError("cannot compute the matrix product between a non-square matrix and SymplecticGivens object."))
    n, k, c, s = (y.form).n, y.k, y.c, y.s
    @inbounds for i in Base.OneTo(n)
        x1, x2 = x[i,k], x[i, n+k]
        x1n, x2n = x[i+n,k], x[i+n,n+k]
        x[i,k] = c * x1 - s * x2
        x[i+n,k] = c * x1n - s * x2n
        x[i,n+k] = s * x1 + c * x2
        x[i+n,n+k] = s * x1n + c * x2n
    end
    return x
end
@inline function LinearAlgebra.lmul!(x::SymplecticGivens{F,N,T}, y::AbstractVecOrMat) where {F<:PairForm,N<:Int,T}
    LinearAlgebra.require_one_based_indexing(y)
    size(y, 1) == size(y, 2) || throw(ArgumentError("cannot compute the matrix product between a SymplecticGivens object and a non-square matrix."))
    n, k, c, s = (x.form).n, x.k, x.c, x.s
    @inbounds for i in Base.OneTo(n)
        y1, y2 = y[2k-1,2i-1], y[2k,2i-1]
        y1n, y2n = y[2k-1,2i], y[2k,2i]
        y[2k-1,2i-1] = c * y1 + s * y2
        y[2k-1,2i] = c * y1n + s * y2n
        y[2k,2i-1] = -s * y1 + c * y2
        y[2k,2i] = -s * y1n + c * y2n
    end
    return y
end
@inline function LinearAlgebra.rmul!(x::AbstractVecOrMat, y::SymplecticGivens{F,N,T}) where {F<:PairForm,N<:Int,T}
    LinearAlgebra.require_one_based_indexing(x)
    size(x, 1) == size(x, 2) || throw(ArgumentError("cannot compute the matrix product between a non-square matrix and SymplecticGivens object."))
    n, k, c, s = (y.form).n, y.k, y.c, y.s
    @inbounds for i in Base.OneTo(n)
        x1, x2 = x[2i-1,2k-1], x[2i-1,2k]
        x1n, x2n = x[2i,2k-1], x[2i,2k]
        x[2i-1,2k-1] = c * x1 - s * x2
        x[2i,2k-1] = c * x1n - s * x2n
        x[2i-1,2k] = s * x1 + c * x2
        x[2i,2k] = s * x1n + c * x2n
    end
    return x
end
Base.:(*)(x::SymplecticGivens, y::SymplecticGivens) = x.k == y.k ? SymplecticGivens(x.form, x.k, x.c * y.c - x.s * y.s, x.c * y.s + x.s * y.c) : Symplectic(x.form, Matrix(x) * Matrix(y))
Base.:(*)(x::SymplecticGivens, y::Symplectic) = x.form == y.form ? Symplectic(x.form, x * y.data) : x * y.data
Base.:(*)(x::Symplectic, y::SymplecticGivens) = x.form == y.form ? Symplectic(x.form, x.data * y) : x.data * y
Base.:(/)(x::SymplecticGivens, y::SymplecticGivens) = x.k == y.k ? SymplecticGivens(x.form, x.k, (x.c * y.c + x.s * y.s)/(y.c^2 + y.s^2), (-x.c * y.s + y.c * x.s)/(y.c^2 + y.s^2)) : Symplectic(x.form, Matrix(x) / Matrix(y))
Base.:(/)(x::SymplecticGivens, y::Symplectic) = x.form == y.form ? Symplectic(x.form, x * inv(y.data)) : x * inv(y.data)
Base.:(/)(x::Symplectic, y::SymplecticGivens) = x.form == y.form ? Symplectic(x.form, x.data * inv(y)) : x.data * inv(y)
Base.:(/)(x::SymplecticGivens, y::AbstractMatrix) = x * inv(y)
Base.:(/)(x::AbstractMatrix, y::SymplecticGivens) = x * inv(y)
Base.:(\)(x::SymplecticGivens, y::SymplecticGivens) = x.k == y.k ? SymplecticGivens(x.form, x.k, (x.c * y.c + x.s * y.s)/(x.c^2 + x.s^2), (x.c * y.s - y.c * x.s)/(y.c^2 + y.s^2)) : Symplectic(x.form, Matrix(x) \ Matrix(y))
Base.:(\)(x::SymplecticGivens, y::Symplectic) = x.form == y.form ? Symplectic(x.form, inv(x) * y.data) : inv(x) * y.data
Base.:(\)(x::Symplectic, y::SymplecticGivens) = x.form == y.form ? Symplectic(x.form, inv(x.data) * y) : inv(x.data) * y
Base.:(\)(x::SymplecticGivens, y::AbstractMatrix) = inv(x) * y
Base.:(\)(x::AbstractMatrix, y::SymplecticGivens) = inv(x) * y