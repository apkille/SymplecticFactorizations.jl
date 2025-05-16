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