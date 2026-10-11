# Chapter 2: The Detailed Outline of Results (printed pp. 12-21, PNG 14-23)

## Core Idea
A complete statement of the model and its analytical results before the derivations: the force-pulse representation, the Oldroyd upper layer, Blake's lower layer, the mean-velocity formula (2.8)/(4.26), the zero-flux no-penetration result, and the Newtonian-mucus convolution comparison.

## Frameworks Introduced
- **The full model in one page**: coordinates x₁ (in-phase beat direction), x₂ (beat, transport, and metachronal direction), x₃ (vertical; impulses at x₃ = 0, interface at −h₀, free surface at H − h₀).
- **Metachronal kinematics**: Λ = cT; cilia one wavelength apart beat in phase; Λ/l cilia per wavelength.
- **Force-pulse representation (2.1), restated as (2.5)**:
  $$\underline{F}=\underline{F}_0\,\delta(x_3)\sum_{n=0}^{\Lambda/l-1}\sum_{m,s=-\infty}^{\infty}\delta(x_2-nl+s\Lambda)\,\delta\!\left(t-\frac{nl}{c}+mT\right)$$
  F₀ is the mean force per unit volume of one cilium's effective stroke; l is the row spacing.
- **Oldroyd upper layer (2.2) and its linearization (2.3)**: $(1+\lambda_1\delta/\delta t)\tau=2\mu_2(1+\lambda_2\delta/\delta t)S$; λ₁ > λ₂; the linearized form uses D/Dt = d/dt + v·∇.
- **Momentum (2.4)**: div τ − ∇p + F = ρ₂ Dv/Dt with continuity div v = 0; stress tensor has only τ₂₃ = τ₃₂ (proved in Appendix A).
- **Lower layer**: Blake's average shear flow plus a parabolic profile from a constant pressure gradient dπ/dx₂; effective depth h̄ = (0.4-0.5)h.
- **Free surface (2.6)**: τ₂₃ = 0 at x₃ = H − h₀.
- **Prescribed mean velocity**: (1/TH)∫∫ v₂ dx₃ dt = V is a constraint, not a boundary condition.
- **Double Fourier expansion** in x₂ (period Λ) and t (period T); the delta-sum force becomes the exponential series (2.7) via Lighthill [10].

## Key Results (equations)
- **Mean (00) velocity profile** (piecewise in x₃): linear slope F₀Λ/(lμ₂) inside the penetration band −h₀ ≤ x₃ ≤ 0; uniform plug above, plus the pressure-gradient parabola.
- **Surface value**: $v_{00}(H-h_0)=F_0\frac{\Lambda}{l}\left(\frac{h}{\mu_1}+\frac{h_0}{\mu_2}\right)+\frac{d\pi}{dx_2}\left(\frac{Hh}{\mu_1}+\frac{H^2}{2\mu_2}\right)$
- **Mean propagation velocity (2.8)**:
  $$V=F_0\frac{\Lambda}{l}\left(\frac{h}{\mu_1}+\frac{h_0}{\mu_2}\right)-F_0\frac{\Lambda}{l}\frac{1}{\mu_2}\frac{h_0^2}{2H}+\frac{d\pi}{dx_2}\left(\frac{Hh}{\mu_1}+\frac{H^2}{3\mu_2}\right)$$
  V is driven mainly by h/μ₁; the h₀ corrections are small; independent of λ₁, λ₂, t₀.
- **Harmonic coefficients (2.9)**: (v₂)_{m0} expressed with cos(b_m(x₃+h₀)) and sin terms; only harmonics m = kΛ/l survive; dispersion relation $b_m^2=-\frac{2\pi i m\rho_2}{\mu_2T}\frac{T+2\pi i m\lambda_1}{T+2\pi i m\lambda_2}$.
- **Stress coefficients**: mean τ₂₃ involves a Heaviside-step jump at x₃ = 0; spatial harmonics (τ₂₃)_{ms} are non-zero only when m + s = kΛ/l: "there is a sum of impulses as a component of stress".
- **No penetration (2.10)-(2.12)**: with the force applied inside the Newtonian layer, the F₀ term's depth integral is zero: no net flow of the upper layer from ciliary action.
- **Newtonian mucus (end of chapter)**: averaging is invalid (λ₁ = λ₂ contradicts λ₁ ≫ t₀ > λ₂); instead solve the unit-impulse fundamental solution, set λ₁ = λ₂, and convolve with the actual force pulse $F(t)=F_t\cos(\pi t/t_0)$ over [−t₀/2, t₀/2]. Resulting coefficients $(v_F)_m=\frac{(t_0/T)\,v_m\cos(\pi m t_0/T)}{1-4m^2(t_0/T)^2}$: velocity proportional to penetration duration, hence "considerably smaller".
- **Wavy surfaces**: results (2.8)-(2.9) remain valid to leading order when interface/free-surface amplitudes are of the order of the relative penetration depth.

## Key Concepts
- **F₀**: force per unit volume of one cilium's effective stroke (per unit area after x₃ integration).
- **Heaviside step u**: used to restrict coefficient functions to the mucus window.
- **Selection rules**: time harmonics only at m = kΛ/l; spatial harmonics only at m + s = kΛ/l.

## Anti-patterns
- Treating the V condition as a boundary condition (it is not needed).
- Applying the impulse/averaging model to Newtonian mucus: the required time ordering collapses.
- Reading the Newtonian-looking mean-stress equation as license to replace mucus by a Newtonian fluid: "this, however, is a mistake".

## Key Takeaways
1. The chapter is self-contained: all headline formulas appear here before derivation.
2. Transport lives on the serous layer: h/μ₁ dominates V.
3. Zero flux without penetration; proportional-to-t₀/T (vanishing) velocity for Newtonian mucus with penetration.
4. The comparison statements (viscoelastic ≫ Newtonian; penetration necessary) are stated as conclusions in Ch. 2 and only proved in Ch. 4.

## Connects To
- **Ch 4a**: derivation of every formula listed here.
- **Ch 4b**: the wavy-surface and Newtonian-convolution details.
