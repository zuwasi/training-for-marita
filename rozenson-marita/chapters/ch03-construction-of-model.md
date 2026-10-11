# Chapter 3: Construction of Model (printed pp. 22-29, PNG 24-31)

## Core Idea
Every idealization is stated and justified: planar beat, shallow short penetration, Oldroyd mucus, impulse force at depth h₀, time-averaging over t₀, and Blake's solution for the lower layer. The model's validity rests on the time-scale ordering λ₂ ≪ t₀ ≪ λ₁ with t₀ ≪ T and λ₁ "large but not too large".

## Frameworks Introduced
- **Geometry (Fig. 1, PNG 24)**: two-layer plane model; cilia rooted in the Newtonian layer (from x₃ = −h − h₀ to −h₀), penetrating the mucus at average depth h₀; F₀ applied at x₃ = 0 in the middle of the active stroke.
- **Small parameters (3.1)-(3.2)**: t₀/T ≪ 1 (penetration much shorter than the beat period); h_p/h ≪ 1 and h_p/H ≪ 1 (penetration much shallower than both layers). Justification: "everyone is agreed that if it does occur, it is of very short duration and very shallow".
- **Full Oldroyd equation (3.3)** with the convected derivative; linearized to (3.4) using D/Dt by neglecting products of small terms. Oldroyd is "flow-independent, and therefore implies no preliminary assumptions concerning the flow"; it generalizes a suspension of elastic particles in a Newtonian fluid.
- **Relaxation-time ordering**: λ₂ ≪ T (gradients die out after the effective stroke); λ₁ large vs T (no transients within a period; stress nearly constant); λ₁ "not too large" (else stresses from successive cilia accumulate and rupture the mucus; Litt [14]: "cilia prefer to transport a secretion of some, but not too much, elasticity").
- **Impulse idealization (3.5)**: F₀ = (0, F₀, 0) applied as δ-pulses at x₃ = 0, spaced l along x₂ with metachronal lag nl/c: valid because t₀ ≪ λ₁ makes two instants within t₀ unresolvable.
- **Time averaging**: g* = (1/t₀)∫₀^{t₀} g dt; optionally weighted by the beat-force distribution (same result for symmetric beats; may improve very asymmetric ones, footnote 4).
- **Generalized functions**: averaging reduces smoothness from C² (a.e.) to C⁰; integration by parts over a period is done in the sense of generalized functions (footnote 6).
- **Layer equations (3.6)-(3.9)**: N-N layer: Oldroyd + momentum + continuity; N layer: τ = 2μ₁S + momentum + continuity, solved by Blake's array solution.

## Key Concepts
- **h_p vs h₀**: h_p is the maximal penetration extent; h₀ the average depth used in the model.
- **"Moment of beat"** (footnote 5): the small vicinity of the beat instant where functions lose smoothness and behave like δ-functions.
- **Effective stroke duration**: about (0.1-0.2)T; planar in the x₂-x₃ plane (justified by Sugino & Naito [23] and Wasserman [24], against the usual 3-D assumption).
- **Elasticity as memory**: the larger the delay λ₁ − λ₂, the more the material "remembers" effort until the next beat adds to it; this is what transfers energy from the short contact to the long transport.

## Assumptions (with justification)
1. Infinite plane double layer, infinite 2-D ciliary array: good for parameter significance and mutual dependence.
2. Quasi-steady state, constant mean velocity V along +x₂.
3. Identical periodic beats, time lag only along x₂ (experimental/computational data).
4. Planar effective stroke (Sugino & Naito [23], Wasserman [24]).
5. Small penetration (3.1)-(3.2) (experimental consensus).
6. Velocity gradients small away from penetration; large gradients exist only near the tip during penetration and cannot be neglected.
7. Products of small terms neglected in (3.3) → (3.4).
8. Force applied at average depth h₀ and mean x₂ position between entry and exit.
9. Lower layer by Blake's solution; boundary conditions identical except at the beat instant, so reuse is legitimate (penetration is short and shallow).
10. Boundary conditions: continuity of velocity and stress at the interface; no-slip at the bottom; stress-free surface; flat surfaces first, waviness reviewed afterwards.

## Anti-patterns
- Do not assume the time-averaging transfers to Newtonian mucus: λ₁ = λ₂ breaks λ₁ ≫ t₀ > λ₂.
- Do not let λ₁ exceed "not too large": mucus rupture.
- Do not ignore the large tip-region gradients during penetration.
- Do not differentiate through the non-smooth stress field (see ch04a).

## Key Takeaways
1. Every subsequent formula is conditional on the ordering λ₂ ≪ t₀ ≪ λ₁, t₀ ≪ T.
2. The impulse postulate is an energy-matching claim: same energy input as the distributed force.
3. Flat surfaces are a first step; Ch. 4.3 shows the flat solution survives O(h₀/H) waviness.
4. Averaging converts the problem into generalized functions; smoothness must be checked a posteriori.

## Connects To
- **Ch 4a**: executes the program for the penetration case.
- **Appendices**: stress-tensor form and vanishing mean corrections.
