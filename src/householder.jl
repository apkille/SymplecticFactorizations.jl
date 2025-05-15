struct SymplecticHouseholder{F<:SymplecticForm,N<:Int,T}
    form::F
    k::N
    P::T
    SymplecticHouseholder(form::F, k::N, P::T) where {F<:SymplecticForm,N<:Int,T} = new{F,N,T}(form, k, P)
end
SymplecticHouseholder(x::SymplecticHouseholder) = x

function householder(form::SymplecticForm, k::Int, v::V) where {V}
    k < form.n || throw(ArgumentError("k must be less than n."))
    length(v) == (form.n - k + 1) || throw(ArgumentError("the length of v must be equal to n - k + 1"))
    vt = transpose(v)
    P = I - (2 / (vt * v)) * v * vt
    return SymplecticHouseholder(form, k, P)
end

function Base.Matrix(x::SymplecticHouseholder{F,N,T}) where {F<:BlockForm,N<:Int,T}
    n, k, P = (x.form).n, x.k, x.P
    M = Matrix{eltype(P)}(I, 2n, 2n)
    copyto!(@view(M[k:n,k:n]), P)
    copyto!(@view(M[n+k:2n,n+k:2n]), P)
    return M
end
function Base.Matrix(x::SymplecticHouseholder{F,N,T}) where {F<:PairForm,N<:Int,T}
    n, k, P = (x.form).n, x.k, x.P
    M = Matrix{eltype(P)}(I, 2n, 2n)
    @inbounds for i in Base.OneTo(n-k+1)
        @inbounds for j in Base.OneTo(n-k+1)
            M[2(i+k)-3, 2(j+k)-3] = P[i,j]
            M[2(i+k)-2, 2(j+k)-2] = P[i,j]
        end
    end
    return M
end
Base.Array(x::SymplecticHouseholder) = Matrix(x)
function Base.AbstractMatrix{T1}(x::SymplecticHouseholder{F,N,T2}) where {F<:BlockForm,N<:Int,T1,T2}
    n, k, P = (x.form).n, x.k, x.P
    M = Matrix{eltype(T1)}(I, 2n, 2n)
    copyto!(@view(M[k:n,k:n]), P)
    copyto!(@view(M[n+k:2n,n+k:2n]), P)
    return M
end
function Base.AbstractMatrix{T1}(x::SymplecticHouseholder{F,N,T2}) where {F<:PairForm,N<:Int,T1,T2}
    n, k, P = (x.form).n, x.k, x.P
    M = Matrix{eltype(P)}(I, 2n, 2n)
    @inbounds for i in Base.OneTo(n-k+1)
        @inbounds for j in Base.OneTo(n-k+1)
            M[2i+1,2j+1] = P[i,j]
            M[2i+2,2j+2] = P[i,j]
        end
    end
    return M
end

Base.checkbounds(x::SymplecticHouseholder, i::Int, j::Int) = (form = x.form; i <= 2 * form.n && j <= 2 * form.n)
Base.checkbounds(x::SymplecticHouseholder, i::Int) = (form = x.form; i <= (2 * form.n)^2)
Base.isequal(x::SymplecticHouseholder, y::SymplecticHouseholder) = x.form == y.form && x.k == y.k && x.P == y.P
Base.isapprox(x::SymplecticHouseholder, y::SymplecticHouseholder) = x.form == y.form && isapprox(x.k, y.k) && isapprox(x.P, y.P)
Base.size(x::SymplecticHouseholder) = (form = x.form; (2 * form.n, 2 * form.n))
Base.size(x::SymplecticHouseholder, dim::Int) = (form = x.form; 2 * form.n)
Base.length(x::SymplecticHouseholder) = prod(size(x))
Base.axes(x::SymplecticHouseholder) = (form = x.form; (Base.OneTo(2 * form.n), Base.OneTo(2 * form.n)))
Base.eltype(x::SymplecticHouseholder) = eltype(x.P)

Base.@propagate_inbounds function Base.getindex(x::SymplecticHouseholder{F,N,T}, i::Int, j::Int) where {F<:BlockForm,N<:Int,T}
    @boundscheck checkbounds(x, i, j)
    n, k, P = (x.form).n, x.k, x.P
    if i < k || n < i < n+k
        if i == j
            return oneunit(eltype(P))
        else
            return zero(eltype(P))
        end
    elseif (k <= i && k <= j) || (n+k <= i && n+k <= j)
        return P[(i % n)-k+1, (j % n)-k+1]
    else 
        return zero(eltype(P))
    end
end
Base.@propagate_inbounds function Base.getindex(x::SymplecticHouseholder{F,N,T}, i::Int) where {F<:BlockForm,N<:Int,T}
    @boundscheck checkbounds(x, i)
    n = (x.form).n
    row = ((i-1) % 2n) + 1
    col = ((i-1) ÷ 2n) + 1
    return getindex(x, row, col)
end
Base.@propagate_inbounds function Base.getindex(x::SymplecticHouseholder{F,N,T}, i::Int, j::Int) where {F<:PairForm,N<:Int,T}
    @boundscheck checkbounds(x, i, j)
    n, k, P = (x.form).n, x.k, x.P
    if i < 2k-1
        if i == j
            return oneunit(eltype(P))
        else
            return zero(eltype(P))
        end
    elseif 2k-1 <= i && 2k-1 <= j
        if i % 2 == j % 2
            return P[cld(i,2)-2, cld(j,2)-2]
        else
            return zero(eltype(P))
        end
    else 
        return zero(eltype(P))
    end
end
Base.@propagate_inbounds function Base.getindex(x::SymplecticHouseholder{F,N,T}, i::Int) where {F<:PairForm,N<:Int,T}
    @boundscheck checkbounds(x, i)
    n = (x.form).n
    row = ((i-1) % 2n) + 1
    col = ((i-1) / 2n) + 1
    return getindex(x, row, col)
end

LinearAlgebra.adjoint(x::SymplecticHouseholder) = x
LinearAlgebra.inv(x::SymplecticHouseholder) = x
Base.copy(x::SymplecticHouseholder) = SymplecticHouseholder(x.form, copy(x.k), copy(x.P))

@inline function Base.copyto!(dest::AbstractMatrix, src::SymplecticHouseholder{F,N,T}) where {F<:BlockForm,N<:Int,T}
    LinearAlgebra.require_one_based_indexing(dest)
    size(dest, 1) == size(dest, 2) || throw(ArgumentError("cannot copy a SymplecticHouseholder object to a non-square matrix."))
    n, k, P = (src.form).n, src.k, src.P
    tP = eltype(P)
    @inbounds for i in Base.OneTo(k-1)
        @inbounds for j in Base.OneTo(n)
            if i == j
                dest[i,j] = oneunit(tP)
                dest[i+n,j] = zero(tP)
                dest[i,j+n] = zero(tP)
                dest[i+n,j+n] = oneunit(tP)
            else
                dest[i,j] = zero(tP)
                dest[i+n,j] = zero(tP)
                dest[i,j+n] = zero(tP)
                dest[i+n,j+n] = zero(tP)
            end
        end
    end
    copyto!(@view(dest[k:n,k:n]), P)
    copyto!(@view(dest[n+k:2n,n+k:2n]), P)
    return dest
end
#=
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
=#