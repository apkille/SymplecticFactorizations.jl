# SymplecticFactorizations

[![Build Status](https://github.com/apkille/SymplecticFactorizations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/apkille/SymplecticFactorizations.jl/actions/workflows/CI.yml?query=branch%3Amain)

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

To compute the symplectic polar decomposition of `S`, which produces a product of a positive-definite symmetric symplectic matrix `P` and orthogonal symplectic matrix `O`, call `polar`:

```julia
julia> F = polar(S)
Polar{Float64, Matrix{Float64}, Matrix{Float64}}
P factor:
4×4 Matrix{Float64}:
  30.7533  -10.957    -25.8009   10.7173
 -10.957     4.93081    9.40483  -4.22713
 -25.8009    9.40483   21.846    -8.72895
  10.7173   -4.22713   -8.72895   4.87129
O factor:
4×4 Matrix{Float64}:
 -0.151294   0.863879     0.480235  -0.0140828
 -0.259956  -0.503439     0.823959   0.00799922
 -0.480235   0.0140828   -0.151294   0.863879
 -0.823959  -0.00799922  -0.259956  -0.503439
```