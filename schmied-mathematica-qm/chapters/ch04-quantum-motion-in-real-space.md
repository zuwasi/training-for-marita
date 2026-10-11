# Chapter 4: Quantum Motion in Real Space

## Core Idea
Put a particle in a box x ∈ [0, a] (V = ∞ outside), truncate to `nmax` states, and work in two exactly related bases: the **momentum basis** (square-well eigenstates, where T̂ is diagonal) and the **finite-resolution position basis** (grid points, where V̂ and x̂ are approximately diagonal). A type-I discrete sine transform (`FourierDST[..., 1]`) converts vectors and matrices between them in O(n log n). That one trick enables mixed-basis Hamiltonians, **split-step/Trotter propagation**, **imaginary-time propagation** for nonlinear (Gross-Pitaevskii) ground states, multi-particle tensor-product grids, and 3D lattices.

## Frameworks Introduced
- **Units (natural-units recipe)**: choose LengthUnit and MassUnit freely; derive `TimeUnit = MassUnit*LengthUnit^2/ℏ`, `EnergyUnit = ℏ/TimeUnit`, so that ℏ = 1. Convert physical inputs with `Quantity[...]/Unit //UnitConvert //N`.
- **Two bases + DST-I**: φn(x) = √(2/a) sin(nπx/a); grid xj = jΔ, Δ = a/(nmax+1); Xnj = √(2/(nmax+1)) sin(πnj/(nmax+1)) is symmetric, orthogonal, and its own inverse.
- **Mixed basis (mixed representation)**: build T̂ in the momentum basis and V̂ on the grid, then transform one into the other's basis.
- **Split-step propagation**: e^{-iΔt H/ℏ} ≈ e^{τV/2} e^{τT} e^{τV/2} (symmetric Trotter, error O(τ³)), with each factor diagonal in its own basis.
- **Imaginary-time propagation**: replace −i(t−t0)/ℏ → −β, renormalize every step, iterate until a fixed point. This works even when H depends on ψ.
- **Tensor-product grids** for N particles (`KroneckerProduct`, `ArrayFlatten`) and d-dimensional arrays with multi-dimensional `FourierDST`.

## Key Concepts
1. **Finite-resolution position basis** |j⟩ = √Δ Σn φn(xj)|n⟩: orthonormal and local, with ϑj(xj') = δjj'/√Δ. Grid values are therefore ψ(xj) = vj/√Δ and the density is |vj|²/Δ.
2. **Energy cutoff**: truncating at nmax keeps only physics below E_nmax = π²ℏ²nmax²/(2ma²).
3. **Locality approximation**: f(x)ϑj(x) ≈ f(xj)ϑj(x) for smooth f, so V̂ ≈ Σ W(xj)|j⟩⟨j|.
4. **Exact operators in the momentum basis**: ⟨n|x̂|n'⟩ = a/2 on the diagonal, and −8ann'/(π²(n²−n'²)²) when n−n' is odd.
5. **Avoid p̂**: its momentum-basis matrix is dense, and p̂²/(2m) ≠ T̂ at finite nmax. Use T̂ directly.
6. **Avoid finite-difference T̂ and p̂**: they are less accurate and make matrix exponentials hard.
7. **Wigner distribution** W(x, k): real phase-space quasi-probability. Here it is evaluated on grid points up to the Nyquist limit, using a piecewise-constant ψ; the result carries a sinc(kΔ) factor.
8. **GP energy functional**: Ekin + Epot + (κ/2)∫|ψ|⁴. The chemical potential is μ = Ekin + Epot + 2Eint.
9. **Exchange symmetry**: the swap operator Π commutes with H, so its ±1 eigenspaces give bosonic (symmetric) or fermionic (antisymmetric) sectors. Contact interaction g·δ(x1−x2) has no effect on fermions.

## Mental Models
- "DST-I is the basis-change button." Applying it once to a vector gives a vector; applying it as a 2D transform to a matrix gives the operator in the other basis (V = X·U·X).
- Diagonal matrices behave like lists of scalars: exp acts elementwise, so propagation is just `Exp[...]*v`.
- Imaginary time is cooling: e^{−βH} filters out every excited component, leaving the ground state.
- Box size `a` controls how much of the wavefunction fits; `nmax` sets the resolution. Converge both.

## Anti-patterns
- `MatrixExp[-I*HP*Δt/ℏ].v0` for long dynamics. It is slow; use split-step instead.
- Building X with `Table[Sin[...]]` for large nmax. Use `FourierDST`.
- Using `Fourier` (periodic) where the box has hard walls. DST-I builds in zero, antisymmetric boundaries.
- Forgetting `Normalize` after each imaginary-time step. The propagation is non-unitary.
- Plotting `v` as if it were ψ. Divide by √Δ (and by Δ^{3/2} in 3D), and pad zeros at x = 0 and x = a.
- Using Vint(xj1, xj2) directly for singular potentials such as Coulomb. Truncate at |x| < Δ.
- Treating every eigenvector of a ≥3-particle Hamiltonian as physical. Some belong to other permutation irreps.

## Code Examples

```mathematica
(* units with ℏ = 1, then grid and operators *)
LengthUnit = Quantity["Meters"]; MassUnit = Quantity["Kilograms"];
TimeUnit = MassUnit*LengthUnit^2/Quantity["ReducedPlanckConstant"] //UnitConvert;
EnergyUnit = Quantity["ReducedPlanckConstant"]/TimeUnit //UnitConvert;
a = 10; m = ℏ = g = 1; nmax = 12;
Δ = a/(nmax+1); xgrid = Range[nmax]*Δ;
TM = SparseArray[Band[{1,1}] -> Range[nmax]^2*π^2*ℏ^2/(2*m*a^2)];
TP = FourierDST[TM, 1];                          (* kinetic, position basis *)
W[x_] = m*g*x; Wgrid = Map[W, xgrid];
VP = SparseArray[Band[{1,1}] -> Wgrid];
HP = TP + VP;
gsP = -Eigensystem[-N[HP], 1,
   Method -> {"Arnoldi", "Criteria" -> "RealPart", MaxIterations -> 10^6}];
ListLinePlot[Join[{{0,0}}, Transpose[{xgrid, gsP[[2,1]]/Sqrt[Δ]}], {{a,0}}]]
```
What it demonstrates: the mixed-basis Hamiltonian for an arbitrary W(x). For the gravity well it gives E0 ≈ 1.864 against the exact 1.85576 (−AiryAiZero[1]·(1/2)^(1/3)). The `-Eigensystem[-H, 1, ...]` trick returns the lowest eigenpair.

```mathematica
(* exact x̂ in the momentum basis, then DST to the grid *)
xM = SparseArray[{Band[{1,1}] -> a/2,
   {n1_,n2_} /; OddQ[n1-n2] -> -8*a*n1*n2/(π^2*(n1^2-n2^2)^2)}, {nmax,nmax}];
xP = FourierDST[xM, 1];      (* ≈ DiagonalMatrix[xgrid] *)
u = FourierDST[v, 1];         (* position → momentum coefficients; self-inverse *)
```
What it demonstrates: converting operators and vectors between bases; that x̂ is almost diagonal on the grid.

```mathematica
propApprox[Δt_?NumericQ, M_Integer /; M >= 1, v0_ /; VectorQ[v0, NumericQ]] :=
 Module[{τ, Ke, Pe2, propKin, propPot2, prop},
  τ = -I*N[Δt/(M*ℏ)];
  Ke = Exp[τ*Range[nmax]^2*π^2*ℏ^2/(2*m*a^2)];
  propKin[v_] := FourierDST[Ke*FourierDST[v, 1], 1];
  Pe2 = Exp[τ/2*Wgrid];
  propPot2[v_] := Pe2*v;
  prop[v_] := propPot2[propKin[propPot2[v]]];
  Transpose[{Range[0, M]/M*Δt, NestList[prop, v0, M]}]]
(* <x>(t) *)
ListLinePlot[{#[[1]], Abs[#[[2]]]^2 . xgrid} & /@ propApprox[50, 1000, vv]]
```
What it demonstrates: **split-step propagation** with no matrices at all. For a time-dependent potential W(x, t), pass `Wt` and evaluate it at t + Δt/(4M) and t + 3Δt/(4M) inside `FoldList[prop, v0, times]`. For GP dynamics, add `+ κ*Abs[v]^2/Δ` to the potential exponent.

```mathematica
groundstate[g_?NumericQ, δβ_?NumericQ, tolerance_:10^-10] :=   (* a = m = ℏ = 1 *)
 Module[{Ke, propKin, propPot2, prop, v0},
  Ke = Exp[-δβ*Range[nmax]^2*π^2/2] //N;
  propKin[v_] := Normalize[FourierDST[Ke*FourierDST[v,1],1]];
  propPot2[v_] := Normalize[Exp[-δβ/2*(Wgrid + g*(nmax+1)*Abs[v]^2)]*v];
  prop[v_] := propPot2[propKin[propPot2[v]]];
  v0 = Normalize@RandomComplex[{-1-I,1+I}, nmax];
  FixedPoint[prop, v0, SameTest -> (Norm[#1-#2] < tolerance &)]]
```
What it demonstrates: a **nonlinear Schrödinger (Gross-Pitaevskii)** ground state by **imaginary-time propagation**. Energies follow from `Ekin = π^2/2*Range[nmax]^2 . Abs[FourierDST[γ,1]]^2`, `Epot = Wgrid . Abs[γ]^2`, and `Eint = (g/2)(nmax+1) Total[Abs[γ]^4]`. In 3D, use `Outer[W, xgrid, xgrid, xgrid]`, the 3D `FourierDST`, and `nn[v_] := v/Norm[Flatten[v]]`.

```mathematica
(* two particles: tensor product + contact interaction + symmetry sectors *)
id = IdentityMatrix[nmax, SparseArray];
TP = KroneckerProduct[T1P, id] + KroneckerProduct[id, T1P];
overlap4 = FourierDST[Table[KroneckerDelta[n1+n2,n3+n4] + KroneckerDelta[n1+n3,n2+n4]
   + KroneckerDelta[n1+n4,n2+n3] - KroneckerDelta[n1,n2+n3+n4] - KroneckerDelta[n2,n1+n3+n4]
   - KroneckerDelta[n3,n1+n2+n4] - KroneckerDelta[n4,n1+n2+n3],
   {n1,nmax},{n2,nmax},{n3,nmax},{n4,nmax}], 1]/(2*a);
HintP = ArrayFlatten[overlap4];
HP[Ω_, g_] = TP + Ω*VP + g*HintP;
Π = ArrayFlatten[SparseArray[{i_,j_,j_,i_} -> 1, {nmax,nmax,nmax,nmax}]];
symmetry[v_] := Re[Conjugate[v].(Π.v)]
```
What it demonstrates: exact quartic overlaps ∫ϑϑϑϑ computed with a 4D DST-I, and the exchange operator. Projectors onto the bosonic and fermionic subspaces come from `Select[Transpose[Eigensystem[Normal[Π]]], #[[1]] == ±1 &][[All,2]] //Orthogonalize`.

## Reference Tables

| Operator | Momentum basis \|n⟩ | Finite-resolution position basis \|j⟩ |
|---|---|---|
| T̂ | exactly diagonal, π²ℏ²n²/(2ma²) | dense (`FourierDST[TM,1]`) |
| V̂ = W(x̂) | dense (`FourierDST[VP,1]`) | ≈ diagonal, W(xj) |
| x̂ | exact analytic matrix (odd n−n') | ≈ diagonal, xj |
| p̂ | dense; avoid | avoid (finite differences are inaccurate) |
| GP term κ\|ψ\|² | n/a | diagonal, κ\|vj\|²/Δ |

| Quantity | Formula / code |
|---|---|
| Grid | Δ = a/(nmax+1), `xgrid = Range[nmax]*Δ` (3D centered: `a*(Range[nmax]/(nmax+1)-1/2)`) |
| ψ(xj), density | vj/√Δ, \|vj\|²/Δ (3D: Δ^{-3/2}, Δ^{-3}) |
| ⟨x̂⟩ | `Abs[v]^2 . xgrid` |
| Gravity-well units | L = (ℏ²/(m²g))^{1/3}, T = (ℏ/(mg²))^{1/3}, E = (mg²ℏ²)^{1/3} |
| Harmonic units | L = √(ℏ/(mω)), E = ℏω |
| BEC interaction | κ = (N−1)·4πℏ²a_s/m |
| Trotter step | real time τ = −iΔt/(Mℏ); imaginary time τ = −δβ |

## Key Takeaways
- `FourierDST[·, 1]` is self-inverse and switches vectors and matrices between the momentum basis and the position basis.
- Build H as `FourierDST[TM,1] + DiagonalMatrix[Wgrid]`; this generalizes to any smooth potential at the cost of slightly lower accuracy.
- Split-step = half potential kick × (DST, kinetic phase, DST) × half potential kick; use `Nest`, `NestList`, or `FoldList`.
- Imaginary-time propagation with `Normalize` plus `FixedPoint` finds nonlinear ground states, where diagonalization cannot.
- Several particles or dimensions: use `KroneckerProduct` / `ArrayFlatten`, or d-dimensional arrays plus multi-dimensional `FourierDST`. The cost grows as nmax^(N·d).
- Converge results in both `a` and `nmax`; errors fall off polynomially, so extrapolate.

## Connects To
- Ch 1: units with `Quantity`, pure functions, `Nest`/`FoldList`/`FixedPoint`, optional arguments.
- Ch 2: basis truncation, `Eigensystem` with Arnoldi, time evolution, tensor products.
- Ch 3: the spin-system techniques reused here for spatial degrees of freedom.
- Ch 5: combining motion with spin, and direct use of p̂.
