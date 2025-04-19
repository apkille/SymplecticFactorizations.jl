# News

## v0.1.6 - 2025-04-19

- Add symplectic Givens support.

## v0.1.5 - 2025-01-30

- Decrease allocations created in `blochmessiah`.

## v0.1.4 - 2025-01-25

- Dispatch base multiplication methods on `Symplectic` type.
- Implement `blochmessiah` decomposition methods and type.

## v0.1.3 - 2025-01-12

- Add `::Type{Symplectic}` to method arguments to create Symplectic matrix types.
- Add support from LinearAlgebra.jl.

## v0.1.2 - 2025-01-10

- Add `Symplectic` matrix type.
- Bump Julia compat to 1.10.0.

## v0.1.1 - 2025-01-07

- Add docstrings for `williamson` and `polar`.
- **(fix)** `O`-`P`-orderings in `polar`.

## v0.1.0 - 2025-01-03

- First release.
- Implement `williamson` and `polar` decomposition methods and objects.
- Add symplectic form types `PairForm` and `BlockForm`.
- Add symplectic tools such as `issymplectic` and `randsymplectic`.