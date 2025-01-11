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