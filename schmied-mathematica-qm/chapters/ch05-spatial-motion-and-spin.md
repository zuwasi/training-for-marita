# Chapter 5: Combining Spatial Motion and Spin

## Core Idea
Glue chapter 3 (spin operators), chapter 4 (finite-resolution position/momentum grids) and chapter 2 (tensor products, `traceout`) into one Hilbert space: motion ⊗ spin. Every operator becomes a `KroneckerProduct` whose factor *order* equals the fixed order of subsystems. You then read physics out of the joint ground state or evolved state through expectation values, reduced density matrices, spin projectors and position projectors. Three worked models: a 1D particle in a spin-dependent force plus transverse field (5.1), a 2D Rashba particle in a box (5.2), and the Jaynes–Cummings atom+cavity in XP phase space (5.3).

## Frameworks Introduced
- **Separable vs non-separable test** (5.1): if every term is diagonal in spin (e.g. `V0(x)+Vz(x) Ŝz`), split into effective Hamiltonians `Ĥ↑,↓ = T + V0 ± Vz/2` and reuse chapter 4 unchanged. A transverse term `Bx Ŝx` mixes sectors, so the joint system must be solved.
- **Tensor-order convention**: 1D: x ⊗ spin. 2D: x ⊗ y ⊗ spin. JC: XP phase space ⊗ pseudo-spin. Trace out "last 2 dims" with `traceout[ψ, -2]`, the motion with `traceout[ψ, nmax]` (or `nmax^2` in 2D).
- **Read-out toolkit**: expectation `Re[Conjugate[ψ].(A.ψ)]`; spin-projected states `(1/2 ± Ŝz)ψ`; position-projected states `(|j⟩⟨j| ⊗ 1)ψ`; normalize local quantities by the projected trace.

## Key Concepts
1. **Spin–motion entanglement**: for `Ĥ = p²/2m + mω²x²/2 − F x Ŝz` (Bx=0) the ground states are Gaussians displaced by ±δ/2, δ = F/(mω²), degenerate with E = ħω/2 − F²/(8mω²); ⟨x̂⊗Ŝz⟩ = δ/4 > 0.
2. **Tunnel splitting**: small Bx>0 picks the antisymmetric combination; splitting ΔE = Bx·exp(−F²/(4mħω³)), exponentially small for strong force.
3. **Correlation without means**: numerically ⟨x⟩ = 0 and ⟨Sz⟩ = 0, yet ⟨x̂⊗Ŝz⟩ ≈ 0.095.
4. **Spin-specific densities** ρ↑, ρ↓ have trace = probability of that spin (0.5 each), sum to ρ(x,x′).
5. **Local spin density matrices** at grid point j have trace = probability at j; local mean spin = `Tr[ρs.sz]/Tr[ρs]`.
6. **Rashba coupling** κ̂z = p̂x⊗Ŝy − p̂y⊗Ŝx, Ĥ = (p̂x²+p̂y²)/2m + V + γŜz + α(p̂xŜy − p̂yŜx); built in the momentum representation; projectors to position transformed with `FourierDST`.
7. **Same matrix, different slot**: `pM` represents both p̂x and p̂y; only its Kronecker position distinguishes them (Q5.2).
8. **Jaynes–Cummings in quadratures**: X̂ = (â+â†)/√2, P̂ = (â−â†)/(i√2), so Ĥ_JC = ħωa Ŝz + ħωc(P̂²/2+X̂²/2) + √2 ħg(X̂⊗Ŝx − P̂⊗Ŝy). Pseudo-spin: Ŝ+ = |e⟩⟨g|, Ŝ− = |g⟩⟨e|.
9. **Wigner negativity** of the reduced cavity state after evolution signals a nonclassical field.

## Mental Models
- A joint state vector is an nmax×2 array flattened row-major: spin is the fast index. That is why spin goes last in every `KroneckerProduct` and `traceout[ψ,-2]` removes it.
- Projector-then-trace = "condition on a measurement outcome": unnormalized, so the trace is the outcome probability.
- The JC cavity is "motion in amplitude space": X̂ is the field amplitude, P̂ its conjugate; no real-space motion.

## Anti-patterns
- Swapping factor order between terms (e.g. `KroneckerProduct[sz, xop]` in one term, x-first elsewhere): silently wrong Hamiltonian.
- Using `P.P` for P̂² on a finite grid; use the exact kinetic form `P2` (In[423]-style) for accuracy.
- Leaving `KroneckerProduct[P,Sy]` complex-typed: wrap with `Re[...]` (imaginary parts are zero by construction).
- Plotting raw coefficients of ψ instead of computing observables/reduced densities.
- Forgetting to divide local spin expectations by `Tr[ρs]`.
- Treating ⟨x⟩=0, ⟨Sz⟩=0 as "uncorrelated".

## Code Examples

```mathematica
(* 5.1: 1D particle with spin in harmonic trap + spin-dependent force + transverse field *)
a = m = ħ = 1; nmax = 100; Δ = a/(nmax+1);
xgrid = a*(Range[nmax]/(nmax+1) - 1/2);
xop = SparseArray[Band[{1,1}] -> xgrid];
idx = IdentityMatrix[nmax, SparseArray];
ids = IdentityMatrix[2, SparseArray];
{sx,sy,sz} = Table[SparseArray[PauliMatrix[i]/2], {i,3}];
TM = SparseArray[Band[{1,1}] -> Range[nmax]^2*π^2*ħ^2/(2*m*a^2)];
TP = FourierDST[TM, 1];
H[ω_, F_, Bx_] = KroneckerProduct[TP, ids] + m*ω^2/2*KroneckerProduct[xop.xop, ids] -
   F*KroneckerProduct[xop, sz] + Bx*KroneckerProduct[idx, sx];
Clear[gs];
gs[ω_?NumericQ, F_?NumericQ, Bx_?NumericQ] := gs[ω, F, Bx] =
  -Eigensystem[-H[N[ω], N[F], N[Bx]], 1,
    Method -> {"Arnoldi", "Criteria" -> "RealPart", MaxIterations -> 10^6}];
γ = gs[100, 5000, 500][[2, 1]];
Xz = Re[Conjugate[γ].(KroneckerProduct[xop, sz].γ)]   (* 0.0954 *)
```
What it demonstrates: x ⊗ spin ordering, sparse Kronecker assembly, memoized lowest eigenpair via negated Arnoldi, correlation observable.

```mathematica
(* reduced, spin-resolved and position-resolved densities *)
ρx = traceout[γ, -2];  ρs = traceout[γ, nmax];          (* {{0.5,-0.233},{-0.233,0.5}} *)
γup = KroneckerProduct[idx, ids/2 + sz].γ;  ρxup = traceout[γup, -2];
γdn = KroneckerProduct[idx, ids/2 - sz].γ;  ρxdn = traceout[γdn, -2];
{Tr[ρxup], Tr[ρxdn]}                                     (* {0.5, 0.5} *)
γx[j_Integer /; 1 <= j <= nmax] :=
  KroneckerProduct[SparseArray[{j, j} -> 1, {nmax, nmax}], ids].γ;
ρsx[j_Integer /; 1 <= j <= nmax] := traceout[γx[j], nmax];
meansx[j_Integer /; 1 <= j <= nmax] := Tr[ρsx[j].sz]/Tr[ρsx[j]];
ListLinePlot[Transpose[{xgrid, Table[meansx[j], {j, nmax}]}]]
```
What it demonstrates: projectors 1/2 ± Ŝz, partial traces from both sides, conditional (local) spin expectation.

```mathematica
(* 5.2: Rashba in momentum representation; x ⊗ y ⊗ spin *)
nmax = 50; idM = IdentityMatrix[nmax, SparseArray]; idS = IdentityMatrix[2, SparseArray];
(* TM, pM: exact kinetic (In[423]) and momentum (In[444]) matrices in momentum basis *)
HM[γ_, α_] = KroneckerProduct[TM, idM, idS] + KroneckerProduct[idM, TM, idS] +
   γ*KroneckerProduct[idM, idM, sz] +
   α*(KroneckerProduct[pM, idM, sy] - KroneckerProduct[idM, pM, sx]);
xM = FourierDST[SparseArray[Band[{1,1}] -> xgrid], 1];   (* x̂ in momentum basis *)
M[j_] := FourierDST[SparseArray[{j, j} -> 1, {nmax, nmax}], 1];
RsP[i_, j_] := traceout[KroneckerProduct[M[i], M[j], idS].γ, nmax^2];
mzP[i_, j_] := Re[Tr[RsP[i, j].sz]/Tr[RsP[i, j]]]
```
What it demonstrates: three-factor Kronecker products, same `pM` in different slots, basis-transformed position projectors.

```mathematica
(* 5.3: Jaynes–Cummings on an XP grid; phase space ⊗ pseudo-spin *)
ħ = 1; a = 10; nmax = Round[a^2/π]; Δ = a/(nmax+1);
X = SparseArray[Band[{1,1}] -> a*(Range[nmax]/(nmax+1) - 1/2)];
P = FourierDST[SparseArray[{n1_, n2_} /; OddQ[n1-n2] ->
     4*I*n1*n2/(a*(n2^2-n1^2)), {nmax, nmax}], 1];
P2 = FourierDST[SparseArray[Band[{1,1}] -> Range[nmax]^2*π^2/a^2], 1];
idX = IdentityMatrix[nmax, SparseArray];
{Sx,Sy,Sz} = Table[SparseArray[PauliMatrix[i]/2], {i,3}]; idS = IdentityMatrix[2, SparseArray];
Ha = KroneckerProduct[idX, Sz];
Hc = KroneckerProduct[X.X/2 + P2/2, idS];
Hint = Sqrt[2]*(KroneckerProduct[X, Sx] - Re[KroneckerProduct[P, Sy]]);
HP[ωa_, ωc_, g_] = ħ*ωa*Ha + ħ*ωc*Hc + ħ*g*Hint;
(* after propagating ψ (section 2.3.4): *)
ρX = traceout[ψ, -2]; ρS = traceout[ψ, nmax];
{Tr[ρS.Sx], Tr[ρS.Sy], Tr[ρS.Sz]}
WignerDistributionPlot[ρX, {-a/2, a/2}]
```
What it demonstrates: cavity as a quadrature grid (not Fock ladder truncation), `P2` instead of `P.P`, `Re` cleanup, atom/field reduced states and Wigner plot.

## Reference Tables

| Model | Tensor order | Spatial basis | Coupling term |
|---|---|---|---|
| 1D + spin (5.1) | x ⊗ s | position grid, T via `FourierDST` | `-F xop⊗sz + Bx idx⊗sx` |
| 2D Rashba (5.2) | x ⊗ y ⊗ s | momentum basis | `α(pM⊗idM⊗sy − idM⊗pM⊗sx)` |
| JC (5.3) | XP ⊗ s | position grid in X | `√2 g (X⊗Sx − P⊗Sy)` |

| Quantity | Code |
|---|---|
| ⟨A⟩ | `Re[Conjugate[γ].(A.γ)]` |
| spin-up part | `KroneckerProduct[idx, ids/2+sz].γ` |
| motion reduced ρ | `traceout[γ, -2]` |
| spin reduced ρ | `traceout[γ, nmax]` (2D: `nmax^2`) |
| variance of A | `Re[Conjugate[γ].(A.A.γ) - (Conjugate[γ].(A.γ))^2]` |

Weak-coupling JC demo: ωa = ωc = 1, g = 0.1, coherent field at ⟨X⟩ = √2 (1 photon), spin down; by t = 100 the Wigner function shows negative regions.

## Key Takeaways
- Fix one subsystem order and use it in every `KroneckerProduct`, projector and `traceout`.
- Separable spin-diagonal Hamiltonians reduce to chapter 4 problems; any transverse spin term forces a joint solve.
- Interpret ψ via observables, spin projectors (1/2 ± Ŝz) and position projectors; projected traces are probabilities.
- Operators for different coordinates can share one matrix; position in the product distinguishes them.
- Prefer exact kinetic forms (`TM`, `P2`) over squared momentum matrices on finite grids.
- Wigner negativity of the reduced field state flags genuine quantum dynamics.

## Connects To
- Ch 2: `KroneckerProduct`, `traceout` (In[256]/In[257]), time propagation (2.3.4).
- Ch 3: spin-1/2 operators from `PauliMatrix`, Dicke basis, S±.
- Ch 4: position grid, `FourierDST` kinetic/momentum operators (In[423], In[444]), Wigner function (In[481]), harmonic oscillator.
