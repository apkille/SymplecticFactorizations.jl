# SymplecticFactorizations

[![Build status (Github Actions)](https://github.com/apkille/SymplecticFactorizations.jl/workflows/CI/badge.svg)](https://github.com/apkille/SymplecticFactorizations.jl/actions)
[![codecov](https://codecov.io/github/apkille/SymplecticFactorizations.jl/graph/badge.svg?token=JWMOD4FY6P)](https://codecov.io/github/apkille/SymplecticFactorizations.jl)

SymplecticFactorizations is a Julia package for computing special decompositions
of symplectic matrices. 

## Usage

To install SymplecticFactorizations.jl, start Julia and run the following command:

```julia
] add SymplecticFactorizations
```

### `Symplectic` matrix type

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

### Symplectic decompositions

To compute the symplectic polar decomposition of `S`, which produces a product of an orthogonal symplectic matrix `O` and positive-definite symmetric symplectic matrix `P`, call `polar`:

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