# Cheatsheet

## The Two Questions and Their Answers
| Question | Answer | Where |
|---|---|---|
| Must cilia penetrate the mucus? | YES: zero net flux without it; depth h₀ nearly irrelevant ("penetration per se") | ch04a (4.32)-(4.33), ch05 |
| Does non-Newtonian mucus matter? | YES: Newtonian substitution costs ~10⁴ in energy; velocity ∝ t₀/T → 0 | ch05, ch04b |

## Master Formulas
- Mean velocity: V = F₀(Λ/l)(h/μ₁ + h₀/μ₂ − h₀²/(2μ₂H)) + (dπ/dx₂)(Hh/μ₁ + H²/(3μ₂))  (4.26)
- Metachrony: Λ = cT; only harmonics m = kΛ/l (time), m + s = kΛ/l (space-time) are forced.
- Dispersion: b_m² = −(2πimρ₂/μ₂T)·(T + 2πimλ₁)/(T + 2πimλ₂)  (4.27)
- Mucus constitutive law: (1 + λ₁δ/δt)τ = 2μ₂(1 + λ₂δ/δt)S; Newtonian when λ₁ = λ₂.
- Newtonian-mucus coefficient: (v_F)_m = (t₀/T)v_m cos(πmt₀/T)/(1 − 4m²(t₀/T)²); no resonance (limit 1/4).

## Reference Numbers
| Quantity | Value |
|---|---|
| Cilium length (airways) | 5-7 μm |
| Cilium diameter / spacing | 0.2-0.4 μm |
| Metachronal wavelength | 20-40 μm |
| Effective stroke | (0.1-0.2)T |
| Blake effective depth h̄ | (0.4-0.5)h |
| Relative penetration ε = h₀/H | ≲ 0.1 (0.75·10⁻² in the numerics; 1055 terms) |
| Mucus viscosity insensitivity | μ₂ = 50 P ≡ 100 P (same V = 0.6) |
| Serous viscosity effect | μ₁ = 0.01/0.02/0.03 P → V = 0.6/0.3/0.2 |
| Energy penalty, Newtonian mucus | ~10⁴ (sum = (π²+8)/16 at A/T = 100, t₀/T = 0.01) |

## Time-Scale Ordering (the model's backbone)
λ₂ ≪ t₀ ≪ λ₁ and t₀ ≪ T; λ₁ large vs T but "not too large" (rupture, Litt).

## Decision Guide
| Situation | Consequence |
|---|---|
| Cilia reach mucus briefly, mucus viscoelastic | Full transport via (4.26); h₀ nearly irrelevant |
| Cilia never reach mucus | ZERO net ciliary flux |
| Mucus Newtonian (λ₁ = λ₂) | Velocity ∝ t₀/T; ≥10× slower; 10⁴ energy cost |
| Negative pressure gradient (plug) | Transport degrades; positive gradient can compensate weak beat |
| Cilia too short / excess serous fluid | Penetration lost: transport failure |
| Next cilia beat upstream (antiplectic) | Energetically optimal: mucus already thinned |

## Anti-Patterns (do not)
- Replace mucus by a Newtonian fluid of equal viscosity.
- Use the Maxwell model (single relaxation time) for penetration problems.
- Set ε = 0 (singular perturbation trap).
- Average over t₀ when λ₁ = λ₂.
- Treat the mean-velocity V condition as a boundary condition.
- Differentiate through the impulse-carrying τ₂₃ field.
- Choose M and Λ/l with common divisors (aliasing; worst M = Λ/l; best M = Λ/l ± 1; sum past the false period).
- Apply the quasi-steady model to start-up/transients.

## Provenance
Model validated against: μ₂-insensitivity (Figs. 6-7), μ₁-scaling (Figs. 8-9), negative τ₂₃ pulse (Figs. 10-11), elastic-slab plateau (Fig. 6), antiplectic thinning (Fig. 5); prior publication: Liron & Rozenson [13], Rozenson & Liron [18].
