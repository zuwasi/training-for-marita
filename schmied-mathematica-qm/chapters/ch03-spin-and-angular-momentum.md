# Chapter 3: Spin and Angular Momentum

## Core Idea
Angular momentum lives in finite-dimensional Hilbert spaces, so it is the cleanest first target for matrix quantum mechanics. Build memoized sparse spin matrices for any S, embed them into multi-spin spaces with `KroneckerProduct`, write the Hamiltonian in an explicit unit system, then diagonalize. The same toolkit covers one electron in a field, the 87Rb hyperfine structure (Breit-Rabi levels, magic field, RWA), the transverse Ising ring (quantum phase transition, entanglement), and quantum circuits (gates, QFT, phase estimation).

## Frameworks Introduced
- **Dicke basis |S,M⟩**, ordered by *decreasing* M (`sz[S]` diagonal = `Range[S,-S,-1]`).
- **Unit-system declaration**: `MagneticFieldUnit`, `EnergyUnit`, `TimeUnit` as `Quantity`; constants become numbers by division. Natural units: `EnergyUnit = MagneticFieldUnit*Quantity["BohrMagneton"]` makes µB = 1.
- **Tensor-product embedding** 𝟙⊗…⊗A⊗…⊗𝟙 (`op[]`); **interaction picture + RWA**; **circuit = matrix product**.

## Key Concepts
1. **Ladder operators**: S±|S,M⟩ = √(S(S+1)−M(M±1))|S,M±1⟩; Sx = (S₊+S₋)/2, Sy = (S₊−S₋)/(2i).
2. **Electron Zeeman**: H = −µB gₑ S·B, gₑ < 0; E± = ∓½µB gₑ|B|, independent of field direction.
3. **87Rb hyperfine (Breit-Rabi)**: H = A_hfs I·J − µB Bz(gI Iz + gS Sz + gL Lz), I=3/2, J=1/2. Zero field: F=2 at +3A/4, F=1 at −5A/4; low field gF ≈ ∓1/2 (Zeeman shifts); high field: |MI,MJ⟩.
4. **Magic field** B0 ≈ 3.22896 G: |1,−1⟩↔|2,1⟩ splitting is first-order field-insensitive (long qubit coherence); ≈ 16A gI/(3µB gS²).
5. **Transition matrix** T_ij = ⟨i|H1|j⟩ in the H0 eigenbasis: zero entries = no direct coupling (selection rules).
6. **RWA**: Δ = ω − (Ei−Ej)/ħ, Ω = √(|T_ij|²/ħ² + Δ²). Resonant dressed states split by |T_ij|; far detuned: ac Zeeman / light shifts ±|T_ij|²/(4ħΔ).
7. **Transverse Ising ring** H = −(b/2)Σ Sx(k) − Σ Sz(k)Sz(k+1): gap ≈ 0 for |b|<1, (|b|−1)/2 for |b|>1, quantum phase transition at b = ±1.
8. **Correlation** C_kk' = ⟨S(k)·S(k')⟩ − ⟨S(k)⟩·⟨S(k')⟩; S=1/2: −3/4 (singlet) ≤ C ≤ +1/4.
9. **Entropy of entanglement**: von Neumann entropy (log₂) of a reduced density matrix.
10. **Exponential wall**: dimension (2S+1)ᴺ; doubling compute adds about one spin.

## Mental Models
- Normalize eigenvectors *before* `Limit[..., Bz -> 0]` or `Bz -> Infinity`.
- Circuit diagrams read left to right; matrix products right to left.
- QFT = CFT mathematically; speedup only if physical gates cost poly(n). Phase estimation dodges exponential input/tomography.

## Anti-patterns
- Calling `Eigensystem` on exact-number matrices for large systems: use `h[N[b]]` (machine precision).
- `Orthogonalize[evec]` on symbolic eigenvectors: explodes. If vectors are real and already orthogonal, use `#/Sqrt[#.#]& /@ evec`.
- Full `IdentityMatrix[d]` in tensor products: pass `SparseArray` as the second argument.
- Naming the spin count `N` (protected): the book uses `n`.
- Brute-force `NDSolve` of the 6.8 GHz-driven problem: needs `MaxStepSize->10^-5`; use the RWA instead.
- Trusting Arnoldi near degeneracy (b ≈ 0): request more eigenstates (m > 2); judge by sums of overlaps.
- Writing `s[x_] = -x Log[2,x]` without `s[0|0.] = 0`: 0 log 0 gives Indeterminate.

## Code Examples

```mathematica
SpinQ[S_] := IntegerQ[2S] && S>=0
splus[0] = {{0}} //SparseArray;
splus[S_?SpinQ] := splus[S] = SparseArray[Band[{1,2}] ->
   Table[Sqrt[S(S+1)-M(M+1)], {M,S-1,-S,-1}], {2S+1,2S+1}]
sminus[S_?SpinQ] := Transpose[splus[S]]
sx[S_?SpinQ] := sx[S] = (splus[S]+sminus[S])/2
sy[S_?SpinQ] := sy[S] = (splus[S]-sminus[S])/(2I)
sz[S_?SpinQ] := sz[S] = SparseArray[Band[{1,1}]->Range[S,-S,-1], {2S+1,2S+1}]
id[S_?SpinQ] := id[S] = IdentityMatrix[2S+1, SparseArray]
```
What it demonstrates: memoized, banded sparse spin-S operators in the Dicke basis; `SpinQ` guards physical spins.

```mathematica
MagneticFieldUnit = Quantity["Gausses"];
EnergyUnit = Quantity["PlanckConstant"]*Quantity["Megahertz"];
TimeUnit = Quantity["Microseconds"];
µBn = Quantity["BohrMagneton"]/(EnergyUnit/MagneticFieldUnit)  (* 1.3996245 *)
ħn = Quantity["ReducedPlanckConstant"]/(EnergyUnit*TimeUnit)  (* 0.15915494 *)
Ix = KroneckerProduct[sx[3/2], id[1/2], id[0]];  (* likewise Iy, Iz *)
Sz = KroneckerProduct[id[3/2], sz[1/2], id[0]];  (* Lx.. with sx[0] last *)
Jx = Sx + Lx; ...; Fx = Ix + Jx; ...
Hhf = A(Ix.Jx+Iy.Jy+Iz.Jz) - µB*Bz*(gI*Iz+gS*Sz+gL*Lz);
hfc = {µB -> µBn, ħ -> ħn,
  A -> Quantity["PlanckConstant"]*Quantity[3.417341305452145,"GHz"]/EnergyUnit,
  gS -> -2.0023193043622, gL -> -0.99999369, gI -> +0.0009951414};
{eval, evec} = Eigensystem[Hhf] //FullSimplify;
nevec = #/Sqrt[#.#] & /@ evec;
NMinimize[eval[[6]] - eval[[7]] - 2 A /. hfc, Bz]   (* magic field: Bz -> 3.22896 *)
```
What it demonstrates: atomic-physics units, composite-spin operators, symbolic Breit-Rabi diagonalization, keeping parameters symbolic and substituting `hfc` at the end.

```mathematica
op[S_?SpinQ, n_Integer, k_Integer, a_?MatrixQ] /;
   1<=k<=n && Dimensions[a] == {2S+1,2S+1} :=
  KroneckerProduct[IdentityMatrix[(2S+1)^(k-1), SparseArray], a,
                   IdentityMatrix[(2S+1)^(n-k), SparseArray]]
sx[S_?SpinQ, n_Integer, k_Integer] /; 1<=k<=n := op[S, n, k, sx[S]]  (* sy, sz same *)
H[S_?SpinQ, n_Integer/;n>=3, b_] := -b/2*Sum[sx[S,n,k], {k,n}] -
   Sum[sz[S,n,k].sz[S,n,Mod[k+1,n,1]], {k,n}]
gs[b_?NumericQ, m_Integer/;m>=1] := gs[b,m] = -Eigensystem[-h[N[b]], m,
   Method -> {"Arnoldi", "Criteria" -> "RealPart", MaxIterations -> 10^6}] //
   Transpose //Sort //Transpose;
s[0|0.] = 0;  s[x_] = -x*Log[2, x];
EE[S_?SpinQ, ψ_] := Total[s /@ Re[Eigenvalues[traceout[ψ, -Length[ψ]/(2S+1)]]]]
```
What it demonstrates: N-spin embedding, ring Hamiltonian, memoized lowest-m eigenstates (negate to get the lowest energies), sorted ascending, entanglement entropy via `traceout` (Ch 2).

```mathematica
op[n_Integer, k_Integer, a_] /; 1<=k<=n && Dimensions[a]=={2,2} :=
  KroneckerProduct[IdentityMatrix[2^(k-1), SparseArray], a, IdentityMatrix[2^(n-k), SparseArray]]
{id,x,y,z} = Table[SparseArray[PauliMatrix[i]], {i,0,3}];
H[n_Integer, k_Integer] /; 1<=k<=n := op[n, k, (x+z)/Sqrt[2]]
RZ[n_Integer, k_Integer, φ_] /; 1<=k<=n := op[n, k, (1+Exp[I φ])/2*id + (1-Exp[I φ])/2*z]
SWAP[n_Integer, {j_Integer, k_Integer}] /; 1<=j<=n && 1<=k<=n && j!=k :=
  (IdentityMatrix[2^n, SparseArray] + X[n,j].X[n,k] + Y[n,j].Y[n,k] + Z[n,j].Z[n,k])/2
P0 = (id + z)/2 //SparseArray;  P1 = (id - z)/2 //SparseArray;
CTRL[n_Integer, c_ /; VectorQ[c,IntegerQ], A_] /; (Unequal@@c) && Min[c]>=1 &&
   Max[c]<=n && Dimensions[A]=={2^n,2^n} :=
  IdentityMatrix[2^n, SparseArray] + Apply[Dot, op[n,#,P1]&/@c].(A - IdentityMatrix[2^n, SparseArray])
CNOT[n_Integer, j_Integer -> k_Integer] := CTRL[n, {j}, op[n, k, x]]
S = CNOT[2, 1->2] . H[2, 1];  S . {1,0,0,0}   (* {1/Sqrt[2],0,0,1/Sqrt[2]}: Bell state *)
```
What it demonstrates: generic gate construction from projectors; `j->k` is purely syntactic pattern sugar; right-to-left gate order.

```mathematica
QFTblock[n_Integer, i_Integer] /; 1<=i<=n :=
  Apply[Dot, Table[CTRL[n, {j}, RZ[n, i, 2π/2^(j+1-i)]], {j, n, i+1, -1}]] . H[n,i]
QFT[n_Integer] /; n>=1 := Apply[Dot, Table[SWAP[n, {i, n+1-i}], {i, 1, n/2}]] .
  Apply[Dot, Table[QFTblock[n, i], {i, n, 1, -1}]]
u = {1}; U[φ_] = {{Exp[2π I φ]}}; U0 = IdentityMatrix[Length[u], SparseArray];
CTRLU[n_Integer, i_Integer, φ_] /; 1<=i<=n :=
  KroneckerProduct[op[n,i,P0], U0] + KroneckerProduct[op[n,i,P1], U[φ]]
t = 4; ψ0 = Flatten[KroneckerProduct[SparseArray[1->1, 2^t], u]] //Normal;
ψ1 = KroneckerProduct[Apply[Dot, Table[H[t,i], {i,t}]], U0] . ψ0;
ψ2[φ_] = Apply[Dot, Table[CTRLU[t, i, 2^(t-i) φ], {i, t, 1, -1}]] . ψ1;
ψf[φ_] = KroneckerProduct[ConjugateTranspose[QFT[t]], U0] . ψ2[φ];
prob[φ_?NumericQ] := Re[Diagonal[traceout[ψf[N[φ]], -Length[u]]]]
mostprobable[φ_?NumericQ] := (Ordering[prob[φ], -1][[1]] - 1)/2^t   (* 0.2 -> 3/16 *)
```
What it demonstrates: QFT assembly (n H gates, n(n−1)/2 controlled phases, n/2 swaps), verified against 2^(−n/2) e^(2πijk/2ⁿ); phase estimation with inverse QFT = `ConjugateTranspose`.

## Reference Tables

| Constant / unit (book values) | Value |
|---|---|
| µB in MHz/G | 1.3996245 |
| ħ in MHz·µs units | 0.15915494 |
| gₑ = gS | −2.0023193043622 |
| gL | −0.99999369 |
| gI (87Rb) | +0.0009951414 |
| A_hfs (87Rb) | 3.417341305452145 GHz |
| Magic field | 3.22896 G |

| Gate | Matrix (single block) |
|---|---|
| X, Y, Z | σx, σy, σz |
| H | (σx+σz)/√2 |
| R_x,y,z(φ) | (1+e^{iφ})/2·𝟙 + (1−e^{iφ})/2·σ |
| SWAP | (𝟙 + XX + YY + ZZ)/2 |
| √SWAP | (3+i)/4·𝟙 + (1−i)/4·(XX+YY+ZZ) |
| CNOT / CCNOT / CSWAP | `CTRL` with x, x, `SWAP` |

| Ising regime | Ground state | Entanglement |
|---|---|---|
| b → +∞ | \|+x⟩^⊗N | 0 |
| b → −∞ | \|−x⟩^⊗N | 0 |
| 0 < b ≪ 1 | (\|↑…↑⟩+\|↓…↓⟩)/√2 | 1 (cat state) |
| −1 ≪ b < 0 | (\|↑…↑⟩−\|↓…↓⟩)/√2 | 1 |

## Key Takeaways
- One sparse, memoized operator library plus `KroneckerProduct` covers atoms, spin chains and qubits.
- Declare units with `Quantity`; keep symbols, substitute numbers (`/. hfc`) last.
- RWA for driven two-level dynamics; brute-force GHz integration is impractical.
- Numerics: `N[b]`, Arnoldi plus sorting, extra eigenstates near degeneracies.

## Connects To
- Ch 1: `SparseArray`, memoization, `...Q` predicates, `Eigensystem`. Ch 2: tensor products, `traceout`, interaction picture. Ch 4+: motional systems reuse units and diagonalization patterns.
