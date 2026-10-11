# Patterns: Using Mathematica for Quantum Mechanics (Schmied)

## Memoized sparse operator library
**When to use**: any spin-S, qubit, or oscillator calculation that reuses the same operators.
**How**: `sz[S_?SpinQ] := sz[S] = SparseArray[Band[{1,1}] -> Range[S,-S,-1], {2S+1,2S+1}]`; same for `splus`, `sx`, `sy`, `id[S] = IdentityMatrix[2S+1, SparseArray]`.
**Trade-offs**: fast and memory-light; cache grows and goes stale, so `Clear` when definitions change. (Ch 1, 3)

## Tensor-product embedding (`op[]`)
**When to use**: composite systems, N spins, particle + spin, multi-particle grids.
**How**: `op[S,n,k,a] := KroneckerProduct[IdentityMatrix[(2S+1)^(k-1),SparseArray], a, IdentityMatrix[(2S+1)^(n-k),SparseArray]]`. Fix one factor order (e.g. x⊗y⊗spin) and use it everywhere.
**Trade-offs**: general and simple; dimension grows as Π nk; wrong ordering gives no error, only wrong physics. (Ch 2, 3, 5)

## Symbolic first, numbers last
**When to use**: physical Hamiltonians with constants (87Rb hyperfine, Zeeman).
**How**: declare `EnergyUnit`, `MagneticFieldUnit`, `TimeUnit` as `Quantity`; keep A, gS, µB symbolic; diagonalize with `Eigensystem // FullSimplify`; substitute `/. hfc` at the end; probe limits with `Assuming[..., Limit/Series]`.
**Trade-offs**: exact insight (limits, magic field via `NMinimize`); symbolic eigenproblems only scale to small matrices. (Ch 1, 3)

## Lowest eigenpairs of large sparse H
**When to use**: ground and low-lying states of spin chains, grids, spin+motion.
**How**: `gs[b_?NumericQ, m_] := gs[b,m] = -Eigensystem[-h[N[b]], m, Method->{"Arnoldi","Criteria"->"RealPart",MaxIterations->10^6}] // Transpose // Sort // Transpose`.
**Trade-offs**: fast and cached; there is no direct "smallest" option (hence the sign trick); unreliable near degeneracies, so request extra states. (Ch 1, 3, 4, 5)

## Numeric-guarded propagator
**When to use**: time-independent H, many time points.
**How**: `U[Δt_?NumericQ] := MatrixExp[-I*H*N[Δt]/ℏ]`; then `ψt = U[t].ψ0`. If [H(t),H(t')] = 0, exponentiate `Integrate[H[s],{s,t0,t}]`.
**Trade-offs**: avoids symbolic blowup; no closed-form t-dependence; dense `MatrixExp` is slow for big H. (Ch 2)

## Interaction picture + RWA
**When to use**: H = H0 + H1(t) with very different energy scales (GHz drive on hyperfine levels).
**How**: diagonalize H0; form Tij = ⟨i|H1|j⟩; either integrate iℏψ' = H̃(t)ψ with `NDSolve`, or use the closed-form two-level RWA with detuning Δ and Rabi Ω.
**Trade-offs**: large time steps / analytic results; RWA ignores far-detuned terms. (Ch 2, 3)

## Partial trace by reshape, contract, flatten
**When to use**: subsystem states, entanglement, measurement probabilities.
**How**: `rdm[ρ,{dA,dB,dC}]` (traces middle slot) via `ArrayReshape` + `TensorContract[...,{2,5}]` + `Flatten`; pure-state version uses `Transpose[P,{1,3,2}].ConjugateTranspose[P]`; `traceout[ψ, ±d]` removes first/last d dims.
**Trade-offs**: exact and general; full state must fit in memory. (Ch 2, 3, 5)

## Gate algebra from Pauli matrices and projectors
**When to use**: building quantum circuits.
**How**: `H = (x+z)/√2`, `SWAP = (𝟙+XX+YY+ZZ)/2`, `CTRL[n,c,A] = 𝟙 + (Dot @@ (op[n,#,P1]&/@c)).(A−𝟙)`; circuit = matrix product read right to left; inverse QFT = `ConjugateTranspose[QFT[t]]`.
**Trade-offs**: no index bookkeeping; 2ⁿ-sized matrices. (Ch 3)

## Mixed-basis Hamiltonian (DST)
**When to use**: a particle in an arbitrary smooth 1D potential.
**How**: `TM = SparseArray[Band[{1,1}] -> Range[nmax]^2 π² ℏ²/(2 m a²)]`; `TP = FourierDST[TM,1]`; `VP = SparseArray[Band[{1,1}] -> W /@ xgrid]`; diagonalize `TP + VP`.
**Trade-offs**: any W(x), no matrix elements to integrate; potential is only approximately diagonal (converge in nmax and box size a). (Ch 4)

## Split-step / Trotter propagation
**When to use**: real-time dynamics, including time-dependent potentials.
**How**: `propKin[v_] := FourierDST[Ke*FourierDST[v,1],1]`; potential half-steps `Exp[-I Δt/2 W]*v`; iterate with `Nest`/`NestList`/`FoldList`.
**Trade-offs**: O(n log n) per step, no matrices; O(τ³) error per step, so many steps. (Ch 4)

## Imaginary-time ground state
**When to use**: non-linear (Gross-Pitaevskii) or very large systems.
**How**: same propagator with τ = −δβ, `Normalize` after every step, `FixedPoint[prop, v0, SameTest -> (Norm[#1-#2] < tol &)]`.
**Trade-offs**: handles ψ-dependent H; slow when the gap is small; ground state only. (Ch 4)

## Symmetry-sector projection
**When to use**: identical bosons/fermions.
**How**: build Π with `SparseArray[{i_,j_,j_,i_}->1, ...]`, select ±1 eigenvectors, `Orthogonalize`, project `s.H.Transpose[s]`.
**Trade-offs**: smaller physical sectors; with ≥ 3 particles other permutation irreps appear. (Ch 4)

## Project, then trace (spin- or position-resolved readout)
**When to use**: spin-resolved densities, local magnetization, conditional expectations.
**How**: `KroneckerProduct[idx, ids/2 + sz].γ` (spin-up part) or a grid-point projector, then `traceout`; divide by `Tr` for conditional values.
**Trade-offs**: result is unnormalized; correlations lost in each reduced matrix. (Ch 5)
