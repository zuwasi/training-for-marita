# Cheatsheet: Using Mathematica for Quantum Mechanics (Schmied)

## Wolfram Language essentials
| Want | Use |
|---|---|
| Evaluate once vs every call | `f[x_] = ...` vs `f[x_] := ...` |
| Rule evaluated now vs per match | `->` vs `:>`; apply `/.`, repeat `//.` |
| Cache results | `f[x_] := f[x] = ...` (and `Clear[f]` before redefining) |
| Guard arguments | `x_?NumericQ`, `x_Integer`, `x_ /; x > 0`, default `b_:5` |
| Local variables | `Module[{i}, ...]`, `With[{c = ...}, ...]` |
| Fast numerics | `N[...]` inputs; never diagonalize exact large matrices |
| Real variables | `Assuming[Element[x, Reals], ...]`, `$Assumptions` |
| Debug | `Echo[x, "label"]`, `Trace`, `FullForm` |
| Sparse identity | `IdentityMatrix[d, SparseArray]` |
| Sorted eigensystem | `Transpose[Sort[Transpose[Eigensystem[H]]]]` |
| Lowest k states | `-Eigensystem[-H, k, Method->{"Arnoldi","Criteria"->"RealPart"}]` |

## Which time-evolution method?
| Situation | Method |
|---|---|
| H time-independent | `MatrixExp[-I H t/ℏ].ψ0` (guard `t_?NumericQ`) |
| [H(t),H(t')] = 0 | `MatrixExp[-I/ℏ Integrate[H[s],{s,t0,t}]].ψ0` |
| H = H0 + H1(t), scales differ | interaction picture + `NDSolve`, or RWA |
| Real-space particle, long times | split-step Trotter with `FourierDST` |
| Ground state of non-linear H | imaginary-time propagation + `Normalize` |
| Truncated series needed | Magnus (unitary), not Dyson |

## Which basis is an operator diagonal in? (Ch 4)
| Operator | Momentum basis (sin modes) | Position grid |
|---|---|---|
| Kinetic T̂ | exact diagonal n²π²ℏ²/(2ma²) | dense (`FourierDST[TM,1]`) |
| Potential V̂ = W(x) | dense | ≈ diagonal W(xj) |
| Position x̂ | analytic, banded-odd | ≈ diagonal xj |
| Momentum p̂ | dense; avoid, p̂²/2m ≠ T̂ | dense |

## Grid rules (Ch 4)
- `Δ = a/(nmax+1); xgrid = Range[nmax] Δ`; ψ(xj) = vj/√Δ; pad zeros at 0 and a for plots.
- `FourierDST[·,1]` is self-inverse; vectors: one call; matrices: V = X.U.X.
- Converge in both `nmax` (resolution, energy cutoff) and box size `a`.
- Natural units: pick length and mass, set ℏ = 1, derive time and energy.

## Composite-system rules
- First `KroneckerProduct` factor = slowest index: |↑↑⟩, |↑↓⟩, |↓↑⟩, |↓↓⟩.
- Product state: `Flatten[KroneckerProduct[ψ1, ψ2]]`; ket-bra: `KroneckerProduct[v, Conjugate[w]]`.
- `traceout[ψ, -2]` drops last qubit/spin; `traceout[ψ, 2]` drops first; `rdm[ρ,{dA,dB,dC}]` drops middle.
- Expectation: `Re[Conjugate[ψ].(A.ψ)]`; probabilities: `Re[Diagonal[ρ]]`.
- Entropy: `s[0|0.] = 0; s[x_] = -x Log[2,x]`; sum over `Eigenvalues[ρA]`.

## Spin and circuits (Ch 3)
- Spin-S: Dicke order M = S, S−1, …, −S; Sx = (S₊+S₋)/2, Sy = (S₊−S₋)/(2i).
- 87Rb: I = 3/2, J = 1/2, A_hfs = 3.417341 GHz, magic field 3.22896 G, µB = 1.3996245 MHz/G.
- Ising ring: gap ≈ 0 for |b| < 1, (|b|−1)/2 for |b| > 1; close ring with `Mod[k+1, n, 1]`; don't name the count `N`.
- Gates act right to left: `CNOT[2,1->2] . H[2,1]` applies Hadamard first (Bell state).
- QFT matrix 2^(−n/2) e^(2πijk/2ⁿ); phase estimation reads φ ≈ `(Ordering[prob, -1][[1]] - 1)/2^t`.
