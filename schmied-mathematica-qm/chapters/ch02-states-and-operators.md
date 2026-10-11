# Chapter 2: Quantum Mechanics: States and Operators

## Core Idea
A computer understands numbers, not Hilbert spaces. Pick an orthonormal basis {|i>}, and every operator becomes a complex matrix `A_ij = <i|A|j>` and every state a complex vector `psi_i = <i|psi>`. Once you do that, the Schrödinger equation turns into linear algebra: `H.psi = E psi` (eigenproblem, Eq. 2.14) or `i hbar psi'(t) = H(t).psi(t)` (coupled ODEs, Eq. 2.25). Composite systems are built with `KroneckerProduct`, and subsystems are extracted by partial trace.

## Frameworks Introduced
- **Basis sets and representations**: use them for every problem. Matrix `A = Sum_ij A_ij |i><j|`, vector `|psi> = Sum_i psi_i |i>`. These two equations convert between the abstract and numeric pictures.
- **Incomplete basis sets**: use when the Hilbert space is infinite. A truncated basis gives `Sum_i |i><i| = P` (a projector), with `Q = 1 - P`. Your model drops the `PAQ`, `QAP` and `QAQ` terms, so check that they are negligible.
- **Variational ground-state calculations**: a truncated-basis ground-state energy is an upper bound. It decreases monotonically as you add basis functions, so use this convergence to judge whether the basis is complete enough.
- **Diagonalization**: use `Eigenvalues`, `Eigenvectors` or `Eigensystem` on the Hamiltonian matrix.
- **Time-dependent Schrödinger equation / propagator**: `|psi(t)> = U(t0,t)|psi(t0)>`. The Dyson series is not unitary when truncated. The **Magnus expansion** `U = exp[Sum_k Omega_k]` stays unitary when truncated and converges faster.
- **Interaction picture**: use when `H = H0 + H1(t)` with H0 diagonalizable and its energy scales large. The basis is `|i(t)> = e^{-i E_i t/hbar}|i>`, with `H_ij(t) = <i|H1(t)|j> e^{-i(E_j-E_i)t/hbar}`. The slow phases allow large ODE time steps.
- **Basis construction**: for a single degree of freedom, use a basis known from the Hilbert space itself (Dicke states |S,M>) or the eigenstates of a simple H0 below an energy cutoff. For coupled degrees of freedom (`H = Sum_k H^(k) + H_int`), use tensor-product bases `|i1,...,iN> = |i1>⊗...⊗|iN>`. The basis size is `Prod n_k` and grows fast.
- **Reduced density matrices [partial trace]**: `rho_AC = Tr_B rho_ABC`. Do it by reshaping to a tensor, contracting the B indices, then flattening again.

## Key Concepts
- **Completeness relation** `Sum_i |i><i| = 1`: inserting it twice gives the matrix elements.
- **Projector P / complement Q**: the subspace kept by the truncated basis, and the part it neglects.
- **Magnus expansion**: `Omega_1 = -(i/hbar) Int H dt1`; `Omega_2` contains `[H(t1),H(t2)]`.
- **MatrixExp**: the true matrix exponential `Sum A^n/n!`.
- **Tensor-product structure** `V = V(1)⊗...⊗V(N)`. A single-subsystem operator is `1⊗...⊗h(k)⊗...⊗1`.
- **Flattened basis ordering**: the last subsystem index runs fastest: |1,1,1>, |1,1,2>, ..., |dA,dB,dC>.
- **Product state vs general state**: a product state has `psi_{i1..iN} = psi(1)_i1 ... psi(N)_iN`.
- **Pure-state shortcut**: for `rho = |psi><psi|`, compute the reduced density matrix directly from psi and never build the d×d rho.

## Mental Models
- Eq. 2.10 (abstract) and Eq. 2.14 (matrix) hold the same information, but only the second one is computable. Writing a problem in a basis is the act of modelling it.
- A truncated basis acts like a spatial low-pass filter. `P_nmax(x,y) = 2 Sum sin(n pi x) sin(n pi y)` is a smoothed delta function with resolution about 1/nmax.
- Pattern specificity orders definitions. Mathematica tries `U[t_?NumericQ]` before `U[t_]`, so the fast numeric rule wins whenever it applies.
- Partial trace is "forgetting a particle". The trace always removes the middle slot of `{dA,dB,dC}`.

## Anti-patterns
- `U[t_] = MatrixExp[-I*H*t/hbar]` with a symbolic t gives a huge, slow expression. Use `U[t_?NumericQ] := MatrixExp[-I*H*N[t]/hbar]` instead.
- Using element-wise `Exp[matrix]` for a propagator: use `MatrixExp`.
- Truncating the Dyson series: this breaks unitarity.
- Adding subsystem matrices of different sizes: pad with `IdentityMatrix[Dimensions[hk]]` first.
- Building ρ = |ψ><ψ| for a large pure state just to trace it: use the vector `rdm`.
- Guessing the `KroneckerProduct` order. It is {|↑↑>, |↑↓>, |↓↑>, |↓↓>}: first factor slowest (Q2.10).

## Code Examples

```mathematica
(* propagator for time-independent H; numeric-only fast definition *)
{sx, sy, sz} = Table[PauliMatrix[i], {i, 3}];
H = Sin[θ] Cos[φ] sx + Sin[θ] Sin[φ] sy + Cos[θ] sz;
U[Δt_?NumericQ] := MatrixExp[-I*H*N[Δt]/ℏ]
ψt = U[0.5] . ψ0
```
What it demonstrates: Eq. 2.34, `psi(t) = exp[-i(t-t0)H/hbar] psi(t0)`. The delayed, NumericQ-guarded definition avoids symbolic blowup. `?U` shows the specific rule listed before the general one.

```mathematica
(* coupled degrees of freedom: embed subsystem operators, product states *)
H1 = KroneckerProduct[h1, IdentityMatrix[Dimensions[h2]], IdentityMatrix[Dimensions[h3]]];
H2 = KroneckerProduct[IdentityMatrix[Dimensions[h1]], h2, IdentityMatrix[Dimensions[h3]]];
H3 = KroneckerProduct[IdentityMatrix[Dimensions[h1]], IdentityMatrix[Dimensions[h2]], h3];
H = H1 + H2 + H3;
A = KroneckerProduct[a1, a2, a3];            (* a1⊗a2⊗a3 *)
ψ = Flatten[KroneckerProduct[ψ1, ψ2, ψ3]]   (* product state *)
```
What it demonstrates: Eqs. 2.37, 2.41 and 2.45. Matrices of the same size can be summed, even when the subsystem dimensions differ.

```mathematica
(* reduced density matrix by partial trace over the middle subsystem B *)
rdm[ρABC_?MatrixQ, {dA_Integer /; dA >= 1, dB_Integer /; dB >= 1, dC_Integer /; dC >= 1}] /;
  Dimensions[ρABC] == {dA*dB*dC, dA*dB*dC} :=
  Flatten[TensorContract[ArrayReshape[ρABC, {dA,dB,dC,dA,dB,dC}], {2,5}], {{1,2},{3,4}}]
(* pure state: never build ρ *)
rdm[ψABC_?VectorQ, {dA_Integer /; dA >= 1, dB_Integer /; dB >= 1, dC_Integer /; dC >= 1}] /;
  Length[ψABC] == dA*dB*dC :=
  With[{P = ArrayReshape[ψABC, {dA,dB,dC}]},
    Flatten[Transpose[P, {1,3,2}] . ConjugateTranspose[P], {{1,2},{3,4}}]]
traceout[ψ_?VectorQ, d_Integer /; d >= 1] /; Divisible[Length[ψ], d] := rdm[ψ, {1, d, Length[ψ]/d}]
traceout[ψ_?VectorQ, d_Integer /; d <= -1] /; Divisible[Length[ψ], -d] := rdm[ψ, {Length[ψ]/(-d), -d, 1}]
```
What it demonstrates: Eq. 2.48 as reshape, contract, flatten. Bipartite cases use a trivial subsystem of dimension 1: `rdm[ρAB,{dA,dB,1}]` gives ρA and `rdm[ρAB,{1,dA,dB}]` gives ρB. `traceout[ψ, d]` removes the first d dimensions when d > 0 and the last |d| when d < 0; the matrix versions work the same way. For a product state, `traceout[ψ,-2] == KroneckerProduct[ψ1, Conjugate[ψ1]]` (Q2.11).

## Reference Tables

| Hamiltonian situation | Propagation method | Mathematica |
|---|---|---|
| Time-independent H | `U = exp[-i(t-t0)H/hbar]` (2.34) | `MatrixExp[-I H N[Δt]/ℏ]` |
| `[H(t),H(t')] = 0` for all t, t' | `U = exp[-(i/hbar) Int_{t0}^{t} H(s) ds]` (2.32) | `MatrixExp` of the integrated matrix |
| General H(t) | Solve `i hbar psi' = H(t).psi` (2.25), or Magnus | `NDSolve` on the coefficient vector |
| `H0 + small H1(t)`, disparate scales | Interaction picture, `H_ij(t)` from Eq. 2.31 | Same ODE, larger time steps |

| Object | Construction | Joint index order |
|---|---|---|
| Operator on subsystem k | `KroneckerProduct[1,...,hk,...,1]` | first factor slowest |
| Product state | `Flatten[KroneckerProduct[ψ1,ψ2,...]]` | last index fastest |
| ket-bra |a><b| | `KroneckerProduct[a, Conjugate[b]]` | (row, column) |
| Tr_B of ρ_ABC | `ArrayReshape` → `TensorContract[{2,5}]` → `Flatten[{{1,2},{3,4}}]` | removes the middle slot |

## Key Takeaways
1. Representing a problem in a basis turns quantum mechanics into numerical matrices and vectors; Eq. 2.14 is the central equation.
2. Truncated bases are projectors. They are variational (energy is an upper bound), so check convergence in nmax.
3. Propagate with `MatrixExp` (never element-wise `Exp`), and guard the definition with `?NumericQ`.
4. Use the interaction picture when H0 has large absolute energies; then only energy differences set the time step.
5. Build composite systems with `KroneckerProduct` and `IdentityMatrix` padding; dimensions multiply.
6. Get a partial trace as reshape, contract, flatten. For pure states, work from the vector directly.

## Connects To
- Ch 1: `MatrixExp`, `Eigensystem`, `ArrayReshape`, `Flatten`, and pattern specificity (sections 1.10.4–1.10.5).
- Ch 3: spin systems, Dicke basis, coupled spins via `KroneckerProduct`, entanglement (Q3.19–Q3.20).
- Ch 4: truncated bases for continuous motion (square well, harmonic oscillator, Eq. 4.28 truncation pitfall).
- Ch 5: combining spin and motional degrees of freedom.
