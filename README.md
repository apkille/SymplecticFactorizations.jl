# SymplecticFactorizations

[![Build status (Github Actions)](https://github.com/apkille/SymplecticFactorizations.jl/workflows/CI/badge.svg)](https://github.com/apkille/SymplecticFactorizations.jl/actions)
[![codecov](https://codecov.io/github/apkille/SymplecticFactorizations.jl/graph/badge.svg?token=JWMOD4FY6P)](https://codecov.io/github/apkille/SymplecticFactorizations.jl)

SymplecticFactorizations is a package for computing special decompositions
of symplectic matrices. 

## Usage

Compute a random symplectic matrix `S` by specifying the symplectic form type:

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

julia> S = randsymplectic(J)
4×4 Matrix{Float64}:
 -26.6597   10.4934    22.6423   -9.53373
 -23.6514    8.08135   19.9027   -7.81124
  -1.76113   0.256758   1.67699   0.322509
 -23.9402    8.67196   19.8991   -8.98003

julia> issymplectic(J, S, atol=1e-10)
true
```
Similar methods exist for `PairForm(n)`, which corresponds to the symplectic form equal to the direct
sum of `[0 1; -1 0]`.

To compute the symplectic polar decomposition of `S`, which produces a product of an orthogonal symplectic matrix `O` and positive-definite symmetric symplectic matrix `P`, call `polar`:

```julia
julia> F = polar(S)
Polar{Float64, Matrix{Float64}, Matrix{Float64}}
O factor:
4×4 Matrix{Float64}:
 -0.151275   0.863869    0.480258  -0.0140919
 -0.259969  -0.503456    0.823944   0.00798633
 -0.480217   0.0140998  -0.151275   0.863892
 -0.823968  -0.0079728  -0.259969  -0.503417
P factor:
4×4 Matrix{Float64}:
  30.7533  -10.957    -25.8008   10.7173
 -10.957     4.9308     9.40484  -4.22713
 -25.8008    9.40484   21.846    -8.72893
  10.7173   -4.22713   -8.72893   4.87128

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

julia> F = williamson(J, V)
Williamson{Float64, Matrix{Float64}, Vector{Float64}}
S factor:
4×4 Matrix{Float64}:
 -0.675114   0.666843   0.586008  -0.437426
  0.733678   0.731394  -0.305086  -0.143899
  0.686757  -0.647903  -1.3332     1.17837
 -0.24263   -0.213865   0.785635   0.722449
symplectic spectrum:
2-element Vector{Float64}:
 0.05191138727617343
 1.8136996243167325

julia> isapprox(F.S * V * F.S', Diagonal(repeat(F.spectrum, 2)), atol = 1e-10)
true
```