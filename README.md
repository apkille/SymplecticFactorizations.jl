# SymplecticFactorizations

[![Build status (Github Actions)](https://github.com/apkille/SymplecticFactorizations.jl/workflows/CI/badge.svg)](https://github.com/apkille/SymplecticFactorizations.jl/actions)
[![codecov](https://codecov.io/github/apkille/SymplecticFactorizations.jl/graph/badge.svg?token=JWMOD4FY6P)](https://codecov.io/github/apkille/SymplecticFactorizations.jl)

SymplecticFactorizations is a Julia package for performing optimized linear algebra computations and special decompositions of symplectic matrices. A symplectic matrix `S` is a `2n x 2n` matrix
that satisfies the condition `SJSᵀ = J`, where `J` is a `2n x 2n` invertible, skew-symmetric matrix called the symplectic form.

## Usage

To install SymplecticFactorizations.jl, start Julia and run the following command:

```julia
] add SymplecticFactorizations
```

## `Symplectic` matrix type

The `Symplectic` matrix type is simply a wrapper around a matrix and its corresponding symplectic form.
To give a taste, compute a random symplectic matrix `S` by specifying the symplectic form type:

```julia
julia> using SymplecticFactorizations

julia> J = BlockForm(2)
BlockForm(2)

julia> symplecticform(J)
4×4 Matrix{Float64}:
  0.0   0.0  1.0  0.0
  0.0   0.0  0.0  1.0
 -1.0   0.0  0.0  0.0
  0.0  -1.0  0.0  0.0

julia> S = randsymplectic(Symplectic, J)
4×4 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 -0.116535    -0.181523  1.15808    0.00240053
  0.137176    -0.498848  0.55193   -1.22289
 -0.853172    -0.578052  0.452435  -0.348686
  0.00544036   0.618131  0.493416  -0.35965

julia> issymplectic(S, atol=1e-10)
true
```
Similar methods exist for `PairForm(n)`, which corresponds to the symplectic form equal to the direct
sum of `[0 1; -1 0]`. The matrix type `Symplectic` is nice for keeping track of the symplectic basis
and performing optimized computations, however, standard Julia matrix types are supported in this package:

```julia
julia> s = randsymplectic(J)
4×4 Matrix{Float64}:
 -7.69543  -41.5856  13.5169    92.976
 -2.04414   -8.3759   5.14638   18.4387
  5.89402   33.3788  -9.05499  -74.8917
 -6.10599  -30.7716   9.99773   68.9329

julia> issymplectic(J, s, atol=1e-10)
true

julia> issymplectic(PairForm(2), s, atol=1e-10)
false
```
In the last line of code, the symplectic check failed because we defined `s` with respect to `BlockForm` rather than `PairForm`.

## Symplectic decompositions

To compute the symplectic polar decomposition of `S`, which produces a product of an orthosymplectic matrix `O` and positive-definite symmetric symplectic matrix `P`, call `polar`:

```julia
julia> F = polar(S)
Polar{Float64, Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}, Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}}
O factor:
4×4 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
  0.123739  -0.181484  0.923535   0.314381
  0.316785  -0.36622   0.177255  -0.856803
 -0.923535  -0.314381  0.123739  -0.181484
 -0.177255   0.856803  0.316785  -0.36622
P factor:
4×4 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
  0.816005     0.243795  -0.187156  -0.00132403
  0.243795     0.926977  -0.131778   0.248883
 -0.187156    -0.131778   1.37965   -0.371624
 -0.00132403   0.248883  -0.371624   1.24352

julia> isapprox(F.O * F.P, S, atol = 1e-10)
true
```

To compute the Williamson decomposition of a positive definite matrix `V`, which finds a congruence relation between a symplectic matrix `S` and a diagonal matrix containing the symplectic eigenvalues `spectrum`, call `williamson`:

```julia
julia> X = rand(4, 4); V = X' * X
4×4 Matrix{Float64}:
 1.39837   1.01419   0.819762  0.644317
 1.01419   1.0662    0.452062  0.525889
 0.819762  0.452062  1.1312    1.02089
 0.644317  0.525889  1.02089   1.09805

julia> F = williamson(Symplectic, J, V)
Williamson{Float64, Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}, Vector{Float64}}
S factor:
4×4 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 -0.67512    0.66685    0.586     -0.437417
  0.733677   0.731394  -0.305086  -0.143899
  0.686747  -0.647894  -1.33318    1.17835
 -0.24263   -0.213864   0.785636   0.722449
symplectic spectrum:
2-element Vector{Float64}:
 0.05191288411625232
 1.8137043652121851

julia> isapprox(F.S * V * F.S', Diagonal(repeat(F.spectrum, 2)), atol = 1e-10)
true
```

To compute the Bloch-Messiah/Euler decomposition of a symplectic matrix `S`, which finds the product of an orthosymplectic matrix followed by a diagonal matrix followed by another orthosymplectic matrix, call `blochmessiah`:

```julia
julia> F = blochmessiah(S)
BlochMessiah{Float64, Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}, Vector{Float64}}
O factor:
4×4 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 -0.636034   -0.492216   0.331745    0.493082
  0.696422   -0.690058   0.0200433   0.195995
 -0.331745   -0.493081  -0.636034   -0.492216
 -0.0200437  -0.195995   0.696422   -0.690058
values:
2-element Vector{Float64}:
 1.0016241730101116
 1.8152715263879216
Q factor:
4×4 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 0.451846  -0.0524922  -0.511355  -0.729107
 0.210612   0.329129   -0.699996   0.597764
 0.511355   0.729107    0.451845  -0.0524917
 0.699996  -0.597764    0.210612   0.32913

julia> D = Diagonal(vcat(F.values, F.values .^ (-1))) # sort singular values in direct sum form
4×4 Diagonal{Float64, Vector{Float64}}:
 1.00162   ⋅        ⋅         ⋅ 
  ⋅       1.81527   ⋅         ⋅ 
  ⋅        ⋅       0.998378   ⋅ 
  ⋅        ⋅        ⋅        0.550882

julia> F.O * D * F.Q # == S
4×4 Symplectic{BlockForm{Int64}, Float64, Matrix{Float64}}:
 -0.116535    -0.181523  1.15808    0.0024005
  0.137176    -0.498848  0.55193   -1.22289
 -0.853173    -0.578053  0.452434  -0.348686
  0.00544102   0.618132  0.493417  -0.35965
```

## `SymplecticGivens` type

A symplectic Givens transformation `G(k,θ)` is defined by the block-ordering matrix as

```math
\begin{bmatrix}
I_{k-1} &    &         &         &   & \\
        & c  &         &         & s & \\
        &    & I_{n-k} &         &   & \\
        &    &         & I_{k-1} &   & \\
        & -s &         &         & c & \\
        &    &         &         &   & I_{n-k}
\end{bmatrix}
```
where `c`, `s` denote `cos(θ)` and `sin(θ)`, respectively. The `SymplecticGivens` type, which is a sparse representation of a symplectic Givens transformation, can be called with `givens(form, k, θ)`, where `form` can be either an object of type `BlockForm` or `PairForm`. Stdlib array operations are supported.

```julia
julia> G = givens(BlockForm(3), 2, pi/4)
SymplecticGivens{BlockForm{Int64}, Int64, Float64}(BlockForm(3), 2, 0.7071067811865476, 0.7071067811865475)

julia> Matrix(G)
6×6 Matrix{Float64}:
 1.0   0.0       0.0  0.0  0.0       0.0
 0.0   0.707107  0.0  0.0  0.707107  0.0
 0.0   0.0       1.0  0.0  0.0       0.0
 0.0   0.0       0.0  1.0  0.0       0.0
 0.0  -0.707107  0.0  0.0  0.707107  0.0
 0.0   0.0       0.0  0.0  0.0       1.0

julia> Matrix(G) * G
6×6 Matrix{Float64}:
 1.0   0.0          0.0  0.0  0.0          0.0
 0.0   2.22045e-16  0.0  0.0  1.0          0.0
 0.0   0.0          1.0  0.0  0.0          0.0
 0.0   0.0          0.0  1.0  0.0          0.0
 0.0  -1.0          0.0  0.0  2.22045e-16  0.0
 0.0   0.0          0.0  0.0  0.0          1.0

julia> Matrix(G * G)
6×6 Matrix{Float64}:
 1.0   0.0          0.0  0.0  0.0          0.0
 0.0   2.22045e-16  0.0  0.0  1.0          0.0
 0.0   0.0          1.0  0.0  0.0          0.0
 0.0   0.0          0.0  1.0  0.0          0.0
 0.0  -1.0          0.0  0.0  2.22045e-16  0.0
 0.0   0.0          0.0  0.0  0.0          1.0
```

## Quick Benchmarks

### 

<details>
  <summary>Matrix Inverse</summary>
<p>

```julia
julia> S = randsymplectic(Symplectic, BlockForm(100));

julia> @btime inv($S.data);
  650.084 μs (9 allocations: 414.30 KiB)

julia> @btime inv($S);
  279.875 μs (12 allocations: 1.22 MiB)
```

</p>
</details>


<details>
  <summary>Symplectic Decompositions</summary>
<p>

```julia
julia> @btime polar(S) setup=(S=randsymplectic(Symplectic, BlockForm(100)));
  6.826 ms (29 allocations: 3.08 MiB)

julia> @btime williamson(J, V) setup=(X=rand(200,200); V=X'*X; J=BlockForm(100));
  22.587 ms (129 allocations: 6.77 MiB)

julia> @btime blochmessiah(S) setup=(S=randsymplectic(Symplectic, BlockForm(100)));
  11.867 ms (64 allocations: 4.68 MiB)
```
<p>
</details>

<details>
  <summary>Version Information</summary>
<p>

```julia
julia> versioninfo()
Julia Version 1.11.1
Commit 8f5b7ca12ad (2024-10-16 10:53 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: macOS (arm64-apple-darwin22.4.0)
  CPU: 8 × Apple M1 Pro
  WORD_SIZE: 64
  LLVM: libLLVM-16.0.6 (ORCJIT, apple-m1)
Threads: 1 default, 0 interactive, 1 GC (on 6 virtual cores)
```
<p>
</details>