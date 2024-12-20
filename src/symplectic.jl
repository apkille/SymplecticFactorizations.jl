abstract type SymplecticForm{N} end

struct BlockForm{N<:Int} <: SymplecticForm{N}
    n::N
end
struct PairForm{N<:Int} <: SymplecticForm{N}
    n::N
end
function Base.show(io::IO, x::SymplecticForm)
    print(io, "$(nameof(typeof(x)))($(x.n))")
end
function symplecticform(f::BlockForm)
    N = f.n
    Omega = zeros(2*N, 2*N)
    @inbounds for i in 1:N, j in N:2*N
        if isequal(i, j-N)
            Omega[i,j] = 1.0
        end
    end
    @inbounds for i in N:2*N, j in 1:N
        if isequal(i-N,j)
            Omega[i, j] = -1.0
        end
    end
    return Omega
end
function symplecticform(f::PairForm)
    N = f.n
    Omega = zeros(2*N, 2*N)
    @inbounds for i in Base.OneTo(N)
        Omega[2*i-1, 2*i] = 1.0
        Omega[2*i, 2*i-1] = -1.0
    end
    return Omega
end

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