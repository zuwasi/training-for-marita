# Chapter 4 (part 1): The Solution, Sections 4.1-4.2 (printed pp. 30-46, PNG 32-48)

## Core Idea
The penetration problem is solved analytically by Fourier-transforming the windowed mucus layer in x₃ and expanding in double Fourier series in (x₂, t); the mean flow gives the master formula (4.26) for V, and the no-penetration comparison (§4.2) proves zero net ciliary flux.

## Frameworks Introduced
- **Full system (4.1)-(4.6)**: Newtonian layer (τ = 2μ₁S, momentum, continuity) + N-N layer (linearized Oldroyd, momentum with the F-force, continuity); Blake's velocity profile in the lower layer.
- **Boundary conditions**: (4.8) velocity continuity at the interface x₃ = −h₀; (4.9) stress continuity there; (4.10) stress-free surface at x₃ = H − h₀; plus the prescribed mean velocity V (a constraint, not a boundary condition).
- **Scalar reduction**: only τ₂₃ is non-zero (Appendix A); the system reduces to
  - (4.12): $\tau_{23}+\lambda_1(\partial_t\tau_{23}+v_2\partial_{x_2}\tau_{23})=\mu_2\partial_{x_3}v_2+\lambda_2\mu_2\partial^2_{x_3\partial t}v_2$; the non-linear term $g\equiv v_2\partial_{x_2}\tau_{23}$.
  - (4.13)-(4.14): momentum components with the δ-series force.
- **Fourier transform with a Heaviside window**: U(y) = (e^{iyh₀} − e^{−iy(H−h₀)})/iy; windowed transforms T₂₃*U etc.; key identity $(iyF)*U = iy(F*U) - e^{iyh_0}f(-h_0) + e^{-iy(H-h_0)}f(H-h_0)$ turns boundary values into unknowns.
- **Split into three coefficient systems** after double Fourier expansion (4.20)-(4.25): (00) mean; (m0) time harmonics; (ms) space-time harmonics. Selection rules: G_{m0} = 0 (integration by parts with periodicity); V_{ms} = 0 for s ≠ 0.
- **Dispersion relation (4.27)**: $b_m^2=-\frac{2\pi i m\rho_2}{\mu_2T}\cdot\frac{T+2\pi im\lambda_1}{T+2\pi im\lambda_2}$.
- **Solutions (4.28)-(4.29)**: (v₂)_{m0} and (τ₂₃)_{m0} as combinations of cos(b_m(x₃+h₀)), sin(b_mx₃) with Heaviside windows; interface value for m = kΛ/l:
  $$(v_2)_{m0}(-h_0)=\frac{F_0\frac{\Lambda}{l}\cos(b_m(H-h_0))}{\frac{2\pi ik\frac{\Lambda}{l}\rho_2}{b_mT}\sin(b_mH)+\frac{\mu_1}{h}\cos(b_mH)}$$

## Key Results
- **Mean velocity profile** (piecewise): linear in x₃ with slope F₀Λ/(lμ₂) within −h₀ ≤ x₃ ≤ 0; uniform above; pressure parabola throughout.
- **Interface mean velocity**: v₀₀(−h₀) = F₀(Λ/l)(h/μ₁) + (dπ/dx₂)(Hh/μ₁).
- **MASTER RESULT (4.26)**:
  $$V=F_0\frac{\Lambda}{l}\left(\frac{h}{\mu_1}+\frac{h_0}{\mu_2}\right)-F_0\frac{\Lambda}{l}\frac{1}{\mu_2}\frac{h_0^2}{2H}+\frac{d\pi}{dx_2}\left(\frac{Hh}{\mu_1}+\frac{H^2}{3\mu_2}\right)$$
  - V is independent of λ₁, λ₂, t₀.
  - Dominant term h/μ₁: serous-layer hydraulics, not mucus rheology, sets the scale (while elasticity enables transport at all; see ch05).
  - As h₀ → 0⁺: V → F₀(Λ/l)(h/μ₁) + pressure term: penetration "per se" suffices; its depth is nearly irrelevant.
- **(τ₂₃)_{ms} (4.30)**: non-zero only for m + s = kΛ/l; equals −F₀(Λ/l)cos(2πsx₃/Λ)u(x₃): δ-like stress spikes at beat instants, "no relaxation of stress and no re-distribution of the mechanical energy input".
- **§4.2 WITHOUT penetration**: force moves into the Newtonian momentum equation; solution (4.31) with flux integrals (4.32)-(4.33); the F₀ and pressure terms contribute ZERO net flux: ciliary activity alone cannot move the mucus. Agreement with Blake's conclusion for Newtonian double layers, but here it is the viscoelastic case that matters and there penetration changes everything.

## Key Concepts
- **Window transform U(y)**: encodes the finite mucus layer in the x₃ Fourier transform.
- **Selection rules**: only m = kΛ/l (time) and m + s = kΛ/l (space-time) harmonics are forced.
- **Almost-everywhere statements**: ∂τ₂₃/∂x₂ = 0 holds only a.e.; τ₂₃ is not smooth (impulse sums), so crossing derivatives is illegitimate.

## Anti-patterns
- **"This, however, is a mistake"**: do not infer from the Newtonian-looking mean-stress equation that mucus can be replaced by a Newtonian fluid; with Newtonian mucus the averaging itself is illegitimate and the equations differ.
- Do not differentiate (4.12) in x₂ to conclude ∂τ₂₃/∂x₂ = 0 everywhere; the swap needs smoothness that τ₂₃ lacks.
- Do not treat the V constraint as a boundary condition.

## Key Takeaways
1. Three decoupled Fourier systems; the mean system carries the transport result.
2. V is set by F₀Λ/l and h/μ₁, with only O(h₀/μ₂) penetration corrections.
3. Time harmonics survive only at the metachronal harmonics m = kΛ/l.
4. No penetration gives identically zero ciliary flux.

## Connects To
- **Ch 2**: same formulas, there without proof.
- **Ch 4b**: wavy-boundary and Newtonian-mucus continuations.
- **Appendices**: τ₂₃-only proof and vanishing mean corrections.
