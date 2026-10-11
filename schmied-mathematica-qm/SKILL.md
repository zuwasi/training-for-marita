---
name: schmied-mathematica-qm
description: "Knowledge base from \"Using Mathematica for Quantum Mechanics: A Student's Manual\" by Roman Schmied (arXiv:1403.7050). Use when writing Wolfram Language code for quantum mechanics: matrix representations of states and operators, spin and angular momentum, coupled spins and qubit circuits (QFT, phase estimation), real-space particle motion with DST grids, split-step and imaginary-time propagation, Gross-Pitaevskii, Wigner functions, spin-orbit and Jaynes-Cummings models, or when studying or referencing the book."
allowed-tools:
  - Read
  - Grep
argument-hint: [topic, framework name, or chapter number]
---

# Using Mathematica for Quantum Mechanics: A Student's Manual
**Author**: Roman Schmied (University of Basel) | **Source**: https://arxiv.org/abs/1403.7050 (v3, 2019) | **Pages**: 164 | **Chapters**: 5 | **Generated**: 2026-10-11

## How to Use This Skill

- **Without arguments**: use the core frameworks below.
- **With a topic** (e.g. `split-step`, `magic field`, `partial trace`): find it in the Topic Index and read that chapter file before answering.
- **With a chapter** (e.g. `ch04`): load that chapter file.
- **Browse**: ask "what chapters do you have?" to see the index.

For any topic not covered in the Core section, read the relevant chapter file first. Code in the chapter files is reconstructed from the PDF text; verify it in a Mathematica kernel before relying on output values.

---

## Core Frameworks & Mental Models

### 1. The governing idea: quantum mechanics = linear algebra in a chosen basis
Pick an orthonormal basis {|i⟩}; every operator becomes a matrix Aij = ⟨i|A|j⟩ and every state a vector ψi = ⟨i|ψ⟩. Then the time-independent Schrödinger equation is `Eigensystem[H]` and the time-dependent one is iℏψ' = H(t).ψ. **Writing a problem in a basis is the act of modelling it.** Mathematica is used as a *descriptive language*, not a high-performance solver; the same ideas translate to Python or Matlab.

### 2. Incomplete basis sets and the variational check
Infinite Hilbert spaces must be truncated: Σ|i⟩⟨i| = P (a projector), Q = 𝟙 − P is neglected. The truncated ground-state energy is an upper bound that decreases monotonically with basis size. **Always converge**: grow `nmax` (and box size `a` in real space) until results stop changing.

### 3. Basis construction for coupled degrees of freedom
Use a natural basis per degree of freedom (Dicke states |S,M⟩ for spins, square-well modes for motion), then tensor products: |i1,…,iN⟩ = |i1⟩⊗…⊗|iN⟩ via `KroneckerProduct`. Single-subsystem operators are 𝟙⊗…⊗h(k)⊗…⊗𝟙 (the book's `op[]`). First factor is the slowest index. Dimension is Π nk: the **exponential wall** (doubling compute adds about one spin).

### 4. Reduced density matrices (partial trace)
ρA = TrB ρAB. The book's `rdm[ρ,{dA,dB,dC}]` reshapes to a tensor, contracts the middle subsystem (`TensorContract[...,{2,5}]`), and flattens; `traceout[ψ, ±d]` removes the first or last d dimensions. For pure states, never build ρ = |ψ⟩⟨ψ|. Use these for entanglement entropy, measurement probabilities, and spin- or position-resolved readout.

### 5. Choose the time-evolution method by Hamiltonian structure
- Time-independent H: `U[t_?NumericQ] := MatrixExp[-I H N[t]/ℏ]`.
- [H(t),H(t')] = 0: exponentiate ∫H dt.
- H = H0 + H1(t) with very different scales: **interaction picture** (remove fast phases of H0), or the **rotating-wave approximation** with detuning Δ and generalized Rabi frequency Ω = √(|Tij|²/ħ² + Δ²).
- Truncated expansions: prefer **Magnus** (unitary) over **Dyson** (not unitary).
- Real-space particles: **split-step Trotter** e^{τV/2}e^{τT}e^{τV/2} (error O(τ³) per step).
- Ground state of a non-linear H: **imaginary-time propagation** (τ = −δβ, `Normalize` each step, `FixedPoint`).

### 6. Two bases + DST-I for real-space motion
Put the particle in a box [0, a], truncate at `nmax`. The **momentum basis** (sin modes) makes T̂ exactly diagonal; the **finite-resolution position basis** (grid xj = jΔ, Δ = a/(nmax+1)) makes V̂ and x̂ approximately diagonal (**locality approximation**). `FourierDST[·,1]` converts vectors and matrices between them in O(n log n) and is self-inverse. Build the **mixed-basis Hamiltonian** `FourierDST[TM,1] + DiagonalMatrix[W /@ xgrid]`. Grid values are ψ(xj) = vj/√Δ. Avoid p̂ and finite differences: p̂²/2m ≠ T̂ at finite nmax.

### 7. Units: symbolic first, numbers last
Declare a unit system with `Quantity` (e.g. `EnergyUnit = Quantity["PlanckConstant"] Quantity["Megahertz"]`, or natural units with ℏ = 1 or µB = 1), convert constants by division, keep parameters symbolic, and substitute a rule list (`/. hfc`) only at the end. Arithmetic on dimensionless numbers is fast; `Quantity` arithmetic is slow.

### 8. Wolfram Language habits that make the physics work
- `=` evaluates now, `:=` on every call; `->` vs `:>` likewise. Use `:=` whenever the body depends on the argument.
- **Memoize** operator libraries: `sz[S_?SpinQ] := sz[S] = SparseArray[...]`. `Clear` before redefining.
- Guard numeric functions with `?NumericQ` and feed `N[...]`; never diagonalize large exact matrices.
- Use `SparseArray` and `IdentityMatrix[d, SparseArray]` everywhere.
- Lowest eigenpairs: `-Eigensystem[-H, k, Method->{"Arnoldi","Criteria"->"RealPart"}]`, then sort. Request extra states near degeneracies.
- The kernel remembers evaluation order, not notebook order.

### 9. Model systems as templates
- **Spin-1/2 in a field / 87Rb hyperfine** (Breit-Rabi, magic field 3.22896 G): composite spin operators, symbolic diagonalization, limits.
- **Transverse Ising ring**: quantum phase transition at b = ±1, cat-state ground state near b = 0, correlations and entanglement entropy.
- **Quantum circuits**: gates from Pauli matrices and projectors (`CTRL`, `SWAP`, `CNOT`), circuits as right-to-left matrix products, QFT, phase estimation.
- **Gravity well, square-well dynamics, Wigner function**: real-space toolkit.
- **Gross-Pitaevskii, two interacting particles, 3D BEC**: non-linear and multi-particle grids, exchange symmetry Π.
- **Particle with spin, Rashba coupling, Jaynes-Cummings**: motion⊗spin tensor products, spin projectors, phase-space dynamics.

### Anti-patterns to flag
- Element-wise `Exp[matrix]` instead of `MatrixExp`.
- Symbolic `U[t_] = MatrixExp[...]` for large H (expression blowup).
- Guessing `KroneckerProduct` order; summing subsystem matrices of different size without identity padding.
- `Orthogonalize` on symbolic eigenvectors; normalize real orthogonal vectors with `#/Sqrt[#.#]&`.
- Forgetting `Normalize` in imaginary time; plotting v instead of v/√Δ.
- Using periodic `Fourier` for hard-wall boxes (use DST-I).
- Naming a variable `N` (protected); `0 Log[0]` without `s[0|0.] = 0`.

---

## Chapter Index

| # | Title | Key Frameworks |
|---|-------|----------------|
| [ch01](chapters/ch01-wolfram-language.md) | Wolfram language overview | immediate vs delayed assignment, rules, memoization, patterns, sparse matrices, Eigensystem/Arnoldi, units |
| [ch02](chapters/ch02-states-and-operators.md) | Quantum mechanics: states and operators | basis representations, incomplete bases, propagators, Magnus, interaction picture, tensor products, `rdm`/`traceout` |
| [ch03](chapters/ch03-spin-and-angular-momentum.md) | Spin and angular momentum | spin operators, 87Rb hyperfine, magic field, RWA, transverse Ising, quantum gates, QFT, phase estimation |
| [ch04](chapters/ch04-quantum-motion-in-real-space.md) | Quantum motion in real space | momentum/position bases, DST-I, mixed basis, Wigner, split-step, Gross-Pitaevskii, imaginary time, interacting particles, 3D |
| [ch05](chapters/ch05-spatial-motion-and-spin.md) | Combining spatial motion and spin | separable vs non-separable H, spin projectors, Rashba coupling, Jaynes-Cummings |

## Topic Index

- **87Rb hyperfine structure / Breit-Rabi** → ch03
- **Arnoldi / lowest eigenstates** → ch01, ch03, ch04
- **Bose-Einstein condensate (3D)** → ch04
- **Contact interaction / two particles** → ch04
- **Debugging (Trace, Echo, FullForm)** → ch01
- **Density matrix, reduced / partial trace** → ch02, ch03, ch05
- **Dicke states / spin operators** → ch02, ch03
- **DST-I / FourierDST** → ch04, ch05
- **Entanglement entropy** → ch03
- **Exchange symmetry (bosons/fermions)** → ch04
- **Factorial, many definitions** → ch01
- **Gravity well / Airy functions** → ch04
- **Gross-Pitaevskii / non-linear Schrödinger** → ch04
- **Imaginary-time propagation** → ch04
- **Incomplete basis / truncation** → ch02, ch04
- **Interaction picture** → ch02, ch03
- **Ising model (transverse field)** → ch03
- **Jaynes-Cummings model** → ch05
- **KroneckerProduct / tensor products** → ch02, ch03, ch05
- **Magic magnetic field** → ch03
- **Magnus / Dyson expansion** → ch02
- **MatrixExp / propagator** → ch02, ch04
- **Memoization** → ch01, ch03
- **Mixed basis Hamiltonian** → ch04
- **Patterns, rules, replacements** → ch01
- **Quantum circuits / gates (CNOT, SWAP, CTRL)** → ch03
- **Quantum Fourier Transform / phase estimation** → ch03
- **Rabi oscillations / RWA / dressed states** → ch03
- **Rashba coupling** → ch05
- **Spin-dependent potential / Stern-Gerlach-like** → ch05
- **Split-step / Trotter** → ch04
- **Sparse arrays** → ch01, ch03
- **Units / Quantity / natural units** → ch01, ch03, ch04
- **Wigner distribution** → ch04, ch05

## Supporting Files

- [glossary.md](glossary.md): all key terms with definitions
- [patterns.md](patterns.md): reusable techniques with code skeletons and trade-offs
- [cheatsheet.md](cheatsheet.md): quick-reference tables (method choice, basis diagonality, grid and composite-system rules)

The book's 17 attached notebooks (embedded in the PDF, viewable in Adobe Reader) include SpinOperators.nb, Rubidium87.nb, IsingModel.nb, QuantumCircuits.nb, GravityWell.nb, WignerDistribution.nb, GroundState.nb, ContactInteraction.nb, 3DBEC.nb, RashbaCoupling.nb, and JaynesCummingsModel.nb.

---

## Scope & Limits

This skill covers the book's content (Mathematica 11 era). Newer Mathematica versions add a built-in QuantumFramework (`QuantumState`, `QuantumOperator`); for that paclet use the `wolfram-quantum-framework` skill. For notebook generation or paper reproduction use `wolfram-mathematica`. The book focuses on descriptive clarity, not large-scale efficiency.
