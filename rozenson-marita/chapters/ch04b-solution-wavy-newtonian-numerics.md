# Chapter 4 (part 2): Sections 4.3-4.5 (printed pp. 47-62, PNG 49-64)

## Core Idea
Three extensions complete the solution: wavy boundaries are handled by perturbation in ε and leave mean values unchanged; Newtonian mucus is handled by convolution with a finite force pulse and gives velocity ∝ t₀/T; the viscoelastic series is summed numerically under strict anti-aliasing rules.

## Frameworks Introduced
- **Wavy-surface model (§4.3)**: interface and free surface perturbed by O(h₀) waves (Fig. 3). Justification: Sanderson & Sleigh observed the interface "follows" the ciliary tips; amplitude beyond penetration depth would contradict penetration itself; free-surface thinning/thickening averages out over t₀ (the ~10× acceleration at the beat makes v₃ briefly non-negligible through continuity).
- **Perturbation expansion in ε**: ε^(NN) = h₀/H, ε^(N) = h₀/h = (H/h)ε^(NN) (footnote 8: different x₃ scales per layer); dimensionless variables (4.34); expansions v₂ = v₂₀ + εv₂₁..., v₃ = εv₃₀..., τ₂₃ = τ₂₃₀ + ετ₂₃₁...
- **Zero-order problem = the plane penetration solution (4.50)**: wavy geometry does not change the leading-order mean transport.
- **O(ε) systems (4.51)-(4.62)**: full perturbed equations in both layers; Appendix B proves all (m0) mean corrections vanish (v₂₁, τ₂₂₀, τ₃₃₀, τ₂₃₁, v₃₀): the mean velocity is unaffected at O(ε).
- **Newtonian mucus via convolution (§4.4)**: λ₁ = λ₂ invalidates the impulse/averaging scheme; solve the unit δ-pulse problem with λ₁ = λ₂ = 0, then convolve with the cosine pulse F(t) = F_t cos(πt/t₀) over [−t₀/2, t₀/2], F_t = (π/2)F₀:
  $$(v_F)_m=\frac{(t_0/T)\,v_m\cos(\pi m t_0/T)}{1-4m^2(t_0/T)^2}$$
  No resonance: as 2mt₀/T → 1 the numerator vanishes too; l'Hôpital gives the finite limit 1/4.
- **Numerical summation (§4.5)**: time grid t_j = jΔt, Δt = T/M; sum the (4.27) series over the period; stopping criterion: stop at the harmonic whose mean-square norm is below ε·V; plotted quantity v₂(x₃,t)/F₀ in SI units (F₀ per unit surface).

## Key Results
- **Flat solution valid to O(h₀/H) for mean values** despite "quite unrealistic" flat surfaces.
- **Singular perturbation warning**: computing the zero-order term at ε = 0 gives the qualitatively different zero-penetration problem; keep h₀ > 0 and take h₀ → 0⁺ only in the results.
- **Newtonian-mucus transport**: velocity carries the factor t₀/T for ANY force time-shape (sine and cosine both give it); vanishes as t₀ → 0; converges in 4-6 harmonics to a depth-uniform profile; at least one order of magnitude below the viscoelastic case even if penetration lasted the whole recovery stroke.
- **Numerics**: ε = 0.75·10⁻² needed 1055 terms; velocity is x₃-uniform except briefly at the beat; oscillations persist 1-2 time increments, then the mucus behaves as an elastic slab; negative velocities during recovery are truncation artifacts; the exaggerated peak-to-mean velocity ratio reflects the momentary-stroke idealization.

## Anti-patterns
- **Aliasing control (the thesis's practical numerical lesson)**:
  - never choose M and Λ/l with common divisors;
  - worst choice M = Λ/l (the exponential is always 1);
  - best choice M = Λ/l ± 1 (leaves ±2πikj/M in the exponent);
  - even with coprime values, a false period appears when one number greatly exceeds the other: the period "shrinks" to Δt·j₀ with smaller amplitude; sum to (M·A)/T harmonics, past the false period, until amplitude is negligible.
- Do not apply the impulse-averaging model to Newtonian mucus (λ₁ = λ₂).
- Do not conflate penetration time t₀ with the effective-stroke duration (they differ).
- Do not use the quasi-steady model for start-up from rest: transients involve large velocity gradients that break the model.

## Key Takeaways
1. Waviness is a second-order concern: mean transport is O(ε)-insensitive.
2. The Newtonian penalty is structural: the t₀/T factor has no viscoelastic counterpart.
3. Numerical credibility of the Fourier series depends on coprime M vs Λ/l and summing past the false period.
4. The mucus behaves as an elastic slab over most of the period; all the action is at the beat.

## Connects To
- **Ch 4a**: base solution being extended.
- **Ch 5**: interpretation, figures, and the energy comparison.
- **Appendix B**: the vanishing O(ε) means used above.
