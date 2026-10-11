---
name: rozenson-marita
description: "Knowledge base from \"The Non-Newtonian Approach to the Muco-Ciliary Transport\" by Marita Rozenson (PhD thesis, Weizmann Institute of Science, 1985). Use when applying Rozenson's frameworks for muco-ciliary transport, mucus rheology, viscoelastic (Oldroyd) fluid modeling, ciliary penetration, Fourier solution of pulsatile two-layer flows, or when referencing this thesis's concepts and results."
allowed-tools:
  - Read
  - Grep
argument-hint: [topic, framework name, or chapter number]
---

# The Non-Newtonian Approach to the Muco-Ciliary Transport
**Author**: Marita Rozenson | **Pages**: ~84 (printed 73 + 6 figure plates) | **Chapters**: 5 + 2 appendices | **Generated**: 2026-10-11

Advisors: Prof. Nadav Liron (Mathematics, Technion) and Prof. Alexander Silberberg (Polymer Research, Weizmann Institute).

## How to Use This Skill

- **Without arguments** — load the core frameworks below for reference
- **With a topic** — ask about `penetration`, `Oldroyd`, `wavy surface`, `numerics`, or another indexed topic; I find and read the relevant chapter file
- **With a chapter** — ask for `ch03` or `4.5`; I load that specific file
- **Browse** — ask "what chapters do you have?" to see the index

When you ask about a topic not covered in Core Frameworks below, I will read the relevant chapter file before answering.

---

## Core Frameworks & Mental Models

### 1. The Two Open Questions (the thesis's spine)
Use these as the organizing questions for the whole subject:
- **Q1 (penetration)**: Is direct mechanical contact between ciliary tips and the mucus layer necessary for transport? Answer: YES. Without penetration, net mucus flux is exactly zero; "the important fact is the penetration per se, and not its depth".
- **Q2 (rheology)**: Does the non-Newtonian (viscoelastic) structure of mucus matter for transport efficiency? Answer: decisively YES. Replacing mucus by a Newtonian fluid of the same viscosity costs a factor of about 10,000 in energy for the same transport velocity.

### 2. Double-Layer "Serous-Mucous" Model
Think of the epithelium as two immiscible fluid layers over an infinite plane:
- Lower **periciliary (serous) layer**: Newtonian, viscosity μ₁, depth h, engulfs the beating cilia.
- Upper **mucus layer (N-N)**: viscoelastic, zero-shear-rate viscosity μ₂, thickness H, sitting like a "blanket" over the ciliary tips.
The ciliary array is the external energy source for the upper layer. You do not need to take a position on whether the mucus layer is permanent or transient; assume only that it persists long enough for a quasi-steady-state model.

### 3. Oldroyd Constitutive Equation for Mucus
$(1+\lambda_1\,\delta/\delta t)\,\tau = 2\mu_2(1+\lambda_2\,\delta/\delta t)\,S$ with relaxation times λ₁ (stresses) and λ₂ (rates of strain), λ₁ > λ₂; Newtonian behavior is the special case λ₁ = λ₂. Prefer Oldroyd over the Maxwell model here because penetration involves two time scales (T ≫ t₀), which a single relaxation time cannot represent. Elasticity means the mucus "memorizes" stress: extra stresses persist after strain rates relax, so successive cilia add their efforts instead of each pushing alone.

### 4. Force-Pulse Postulate
Represent each ciliary penetration as a delta-function impulse in time (and in x₂, x₃) applied at the average penetration depth h₀, with magnitude F₀ chosen to match the real energy input. Justified by the time-scale ordering t₀ ≪ λ₁: two instants within the penetration window cannot be resolved.

### 5. Time-Averaging Hierarchy
Average variables over the penetration interval t₀ (uniformly, or weighted by the beat-force distribution). Valid ordering: λ₂ ≪ t₀ ≪ λ₁ with t₀ ≪ T. This averaging destroys smoothness (C² becomes C⁰), so work with generalized functions and their Fourier transforms ("the main mathematical tool used in this work"). CAUTION: the averaging is illegitimate for Newtonian mucus (λ₁ = λ₂); there, use the fundamental solution plus convolution with the actual force pulse instead.

### 6. The Master Result: Mean Propagation Velocity
$$V = F_0\frac{\Lambda}{l}\left(\frac{h}{\mu_1} + \frac{h_0}{\mu_2} - \frac{h_0^2}{2\mu_2 H}\right) + \frac{d\pi}{dx_2}\left(\frac{Hh}{\mu_1} + \frac{H^2}{3\mu_2}\right)$$
Key properties: V is independent of λ₁, λ₂ and t₀; the dominant transport term is h/μ₁ (serous layer), so mucus viscosity μ₂ is of minor importance as long as the mucus is elastic; computed results agree with figures: μ₂ = 50 P vs 100 P give identical transport, while μ₁ = 0.01, 0.02, 0.03 P scale V as 0.6, 0.3, 0.2.

### 7. Zero-Flux-Without-Penetration Result
When the force acts only inside the Newtonian layer (cilia never reach the mucus), the ciliary force term integrates to zero net flux ((4.33): the F₀ term's depth integral vanishes identically). Any observed transport must then come from other mechanisms entirely.

### 8. Elastic-Slab Picture
After the beat, velocity gradients die out within about one relaxation time and the mucus moves as a nearly uniform "elastic slab". Successive effective strokes pile up a "bunch" of mucus; slow stress relaxation spreads each cilium's energy input over the whole period instead of dissipating it locally.

### 9. Antiplectic Synchronization Argument (Fig. 5)
During an effective stroke, cilia thin the mucus above themselves. The next cilia to beat should sit BEHIND (upstream of) the ones that just finished: there the mucus is already thinned, so less energy is wasted re-thinning it. This explains why antiplectic metachronism dominates mucociliary epithelia.

### 10. Anti-Patterns (what this thesis explicitly warns against)
- **Substituting Newtonian mucus**: wastes ciliary energy on local displacements; transport drops ≥10×, energy cost rises ~10⁴.
- **Maxwell single-relaxation-time mucus**: cannot handle two time scales; Ross & Corrsin's "no penetration needed" conclusion "should be reviewed".
- **Deep-penetration slender-body assumptions** (Blake): no experimental support; deep penetration of viscoelastic mucus is energetically unreasonable.
- **Setting ε = 0 in the wavy-boundary problem**: a singular perturbation; the zero-penetration problem is qualitatively different from h₀ → 0⁺.
- **Differentiating across non-smooth stress fields**: τ₂₃ contains impulse sums; ∂τ₂₃/∂x₂ = 0 holds only almost everywhere.
- **Numerical aliasing**: never choose M (time steps per period) and Λ/l with common divisors; worst case M = Λ/l; best M = Λ/l ± 1; sum past the false period ((M·A)/T harmonics) until amplitude is negligible.
- **Treating the mean-velocity condition as a boundary condition**: it is not needed to solve the problem.
- **Using the quasi-steady model for start-up/transients**: it describes only the established periodic regime.

---

## Chapter Index

| # | Title (printed pages) | Key Frameworks |
|---|----------------------|----------------|
| [ch01](chapters/ch01-introduction.md) | Introduction (4-11) | beat commonality, model landscape, two open questions |
| [ch02](chapters/ch02-outline-of-results.md) | Detailed Outline of Results (12-21) | full model summary, V formula, no-penetration result |
| [ch03](chapters/ch03-construction-of-model.md) | Construction of Model (22-29) | force-pulse postulate, averaging hierarchy, assumptions |
| [ch04a](chapters/ch04a-solution-penetration.md) | The Solution: with/without penetration (30-46) | Fourier transform + window, three coefficient systems, (4.26), zero flux |
| [ch04b](chapters/ch04b-solution-wavy-newtonian-numerics.md) | The Solution: wavy, Newtonian, numerics (47-62) | perturbation in ε, convolution, aliasing control |
| [ch05](chapters/ch05-discussion.md) | Discussion (63-67) | energy 10⁴ result, antiplectic argument, failure modes |
| [appendices](chapters/appendices.md) | Appendices A-B + References (68-73) | stress-tensor proof, vanishing O(ε) means, bibliography |

PNG page = printed page + 2 (PNG 1-3 are title/dedication pages without numbers).

## Topic Index

- **Antiplectic synchronization** → ch05, ch01
- **Averaging (time, over t₀)** → ch03, ch04a
- **Boundary conditions** → ch03, ch04a
- **Convolution / fundamental solution** → ch04b, ch02
- **Energy comparison (10⁴ factor)** → ch05
- **Fourier series / transform / generalized functions** → ch04a, ch03
- **Metachronal wave (Λ = cT)** → ch01, ch03
- **No penetration → zero flux** → ch04a, ch02
- **Numerical summation, aliasing** → ch04b
- **Oldroyd equation** → ch03, ch02
- **Penetration (h₀, t₀, h_p)** → ch03, ch04a, ch05
- **Pressure gradient dπ/dx₂** → ch04a, ch05
- **Relaxation times (λ₁, λ₂)** → ch03, ch05
- **Serous vs mucus viscosity roles (h/μ₁ dominance)** → ch04a, ch05
- **Slender body / envelope / Blake models (critique)** → ch01
- **Stress tensor form (τ₂₃ only)** → appendices, ch03
- **V formula (mean velocity)** → ch04a, ch02
- **Wavy surfaces (perturbation in ε)** → ch04b

## Supporting Files

- [glossary.md](glossary.md) — all key terms with definitions
- [patterns.md](patterns.md) — all solution techniques and methods
- [cheatsheet.md](cheatsheet.md) — quick reference tables and decision guides

---

## Scope & Limits

This skill covers the thesis content only. The model is quasi-steady (no transients), planar-beat, infinite-plane, with flat surfaces justified to O(h₀/H); the author notes the equations are valid only within the assumed ranges of λ₁ and λ₂ and that sensitivity to parameters cannot be investigated from the final equations alone. For clinical context or later experimental literature, combine with external sources.
