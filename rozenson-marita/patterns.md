# Patterns (Solution Techniques and Methods)

## Force-Pulse Postulate
**When to use**: modeling short, shallow, periodic contacts between an actuator (cilium) and a medium.
**How**: replace each contact by a δ-impulse in time at the average depth and position, with magnitude matched to the real energy input; sum over the array with metachronal lag δ(t − nl/c) and spacing l.
**Trade-offs**: loses intra-penetration detail (unresolvable anyway when t₀ ≪ λ₁); gains a tractable forcing term; the period-edge stress impulses are artifacts of this choice.

## Time-Averaging over the Contact Interval
**When to use**: the contact duration t₀ is much shorter than the flow period T and than the stress relaxation time λ₁.
**How**: replace fields by g* = (1/t₀)∫₀^{t₀} g dt; optionally weight by the beat-force distribution (equivalent for symmetric beats).
**Trade-offs**: valid ordering λ₂ ≪ t₀ ≪ λ₁; illegitimate for Newtonian mucus (λ₁ = λ₂); averaging downgrades smoothness C² → C⁰, forcing generalized functions.

## Double Fourier Series + Three-System Split
**When to use**: periodic forcing in x₂ (period Λ) and t (period T) with a δ-sum source.
**How**: expand every field in exp(2πisx₂/Λ)exp(2πimt/T) (Lighthill [10] representation of the δ-sum); split into (00), (m0), (ms) coefficient systems; use selection rules (only m = kΛ/l forced in time; only m + s = kΛ/l in space-time).
**Trade-offs**: decouples the problem; (ms) system is unwieldy and was left unsolved (not needed for means).

## Fourier Transform with a Heaviside Window
**When to use**: solving a PDE on a finite layer (x₃ ∈ [−h₀, H−h₀]) with boundary data at the edges.
**How**: multiply fields by u(x₃+h₀) − u(x₃−H+h₀); transform in x₃; convert boundary values into the transform via (iyF)*U = iy(F*U) − e^{iyh₀}f(−h₀) + e^{−iy(H−h₀)}f(H−h₀).
**Trade-offs**: boundary values become unknowns solved by the free-surface and interface conditions; requires "well behaved at infinity" fields.

## Fundamental Solution + Convolution (for the linear/Newtonian case)
**When to use**: when the impulse/averaging scheme is invalid (Newtonian mucus, λ₁ = λ₂) but the problem is linear.
**How**: solve the unit δ-pulse problem (F₀ = 1) at λ₁ = λ₂, then convolve with the actual force pulse F(t) = F_t cos(πt/t₀) over [−t₀/2, t₀/2]; coefficients (v_F)_m = (t₀/T)v_m cos(πmt₀/T)/(1 − 4m²(t₀/T)²).
**Trade-offs**: exact for any force time-shape (only the t₀/T factor survives); depends on linearity (the "duality" of averaged-singular vs continuous-exact forcing, footnote 2).

## Boundary Perturbation in ε (wavy surfaces)
**When to use**: boundary amplitude O(h₀) on a layer of depth H (ε = h₀/H ≲ 0.1).
**How**: nondimensionalize per layer (4.34); expand v, τ, p in ε; zero order = the plane solution (4.50); first-order systems (4.51)-(4.62); prove mean corrections vanish (Appendix B identity (v₂₀v₃₀)_{m0} = 0).
**Trade-offs**: flat-solution validity to O(h₀/H); NEVER set ε = 0 (singular perturbation: the zero-penetration problem is qualitatively different from h₀ → 0⁺).

## Separation of Variables for Vanishing Components (Appendix A)
**When to use**: proving that unwanted stress/velocity components vanish.
**How**: from the homogeneous linearized Oldroyd equations, separate τ = T(x₃,t)X(x₂,x₃); both factors are exponentials in x₂ and t; periodicity in x₂ and t forces the amplitude C(x₃) = 0.
**Trade-offs**: relies on periodicity; scan-level sign ambiguities in the printed exponents do not affect the conclusion.

## Series Summation with Aliasing Control (§4.5)
**When to use**: summing the double Fourier series numerically over a beat period.
**How**: grid t_j = jΔt, Δt = T/M; keep M and Λ/l coprime, ideally M = Λ/l ± 1; sum to (M·A)/T harmonics, past the false period; stop when a harmonic's mean-square norm drops below ε·V.
**Trade-offs**: 1055 terms at ε = 0.75·10⁻²; residual negative-velocity lobes are truncation artifacts; worst-case M = Λ/l makes every exponential equal 1.

## Energy Comparison via Velocity-Squared Bounds
**When to use**: comparing transport efficiency of two rheologies under equal energy input.
**How**: compare kinetic energy per unit volume (velocity squared); bound the Newtonian series v_F² ≤ 2(t₀/T)²v²(1 + Σ(1/(1−4k²(t₀/T)²(A/T)²))²); evaluate the sum in closed form ((π²+8)/16 for A/T = 100, t₀/T = 0.01).
**Trade-offs**: gives the headline 10⁴ factor; depends on the δ-force idealization for the time shape.

## Reuse of a Valid Sub-Solution (Blake's lower layer)
**When to use**: one layer's physics is unchanged by the coupling except at instants that measure zero.
**How**: verify the boundary conditions are identical except at the beat instant; import Blake's array solution and parabolic pressure profile with effective depth h̄ = (0.4-0.5)h unmodified.
**Trade-offs**: legitimate only because penetration is short and shallow; the deep-penetration variant (needed for slender-body mucus solutions) is rejected as experimentally unsupported and energetically unreasonable.
