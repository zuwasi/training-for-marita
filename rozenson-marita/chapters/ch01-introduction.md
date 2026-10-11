# Chapter 1: Introduction (printed pp. 4-11, PNG 6-13)

## Core Idea
Mucociliary transport is clinically essential (lung clearance), the epithelium's structure is essentially the same across organs so a unified model is justified, and existing model families (envelope, slender body, cilia sublayer, Maxwell, Newtonian mucus) each fail for the real two-layer, viscoelastic system.

## Frameworks Introduced
- **Unified epithelium approach**: treat the tissue as a cellular layer covered by long flexible rods (cilia) beating in serous fluid, with mucus on top. Use when modeling any ciliated organ's transport.
- **Model landscape critique** (what each prior approach can and cannot do):
  - **Envelope model**: replaces the dense ciliary array by the wavy surface of the ciliary tips. Valid because cilia are dense; loses cilia-scale physics.
  - **Slender body theory / superposition** (Liron & Meyer [11], Liron & Mochon [12]): Stokes flow built from force distributions along cilia. Fails when the mucus constitutive equation is non-linear (superposition is lost).
  - **Cilia sublayer model** (Blake [2]): singular external forces distributed over the whole beating array. Reused unmodified for the lower Newtonian layer.
  - **Maxwell fluid mucus** (Ross & Corrsin [18]): single relaxation time, wavy-surface cilia. Inadequate: penetration involves two time scales (T ≫ t₀); their "penetration unnecessary" conclusion "should be reviewed".
  - **Newtonian high-viscosity mucus** (Blake [2], Blake & Winet [4], Yates et al. [26]): treats mucus as a viscous fluid; misses elasticity entirely.

## Key Concepts
- **Cilium**: nine outer fibre pairs around one inner pair; bends by relative fibre sliding. Length 5-7 μm in mammalian respiratory epithelium; diameter ~0.2-0.4 μm, comparable to inter-cilium spacing.
- **Effective stroke**: fast, essentially rigid movement of a straightened cilium in the beat direction; the "pushing" phase; also called the active stroke.
- **Recovery stroke**: the return movement; the cilium stays completely submerged in serous fluid.
- **Metachronal wave**: phase-lag wave over the array; wavelength Λ = distance between nearest in-phase cilia; 20-40 μm; Λ = cT.
- **Symplectic / antiplectic**: wave propagating with / against the beat direction. Antiplectic is the most common in mucociliary epithelia; the thesis partly explains why.
- **Serous-mucous structure**: two-layer medium; lower Newtonian layer engulfs the cilia, upper layer has non-linear rheology that must be taken into account.
- **Penetration**: ciliary tips enter the mucus briefly during the effective stroke; experiment shows no (or insufficient) propulsion without it (Cheung & Jahn [6]).

## Mental Models
- Think of the mucus as a "blanket" carried by tips that briefly poke into it; the debate over whether the layer is permanent (Yates et al. [26]) or a time-limited disturbance response (Spungin & Silberberg [22]) need not be resolved: assume only quasi-steady persistence.
- Think of the ciliary array as a conglomerate of long thin bodies, close-packed, submerged in Newtonian liquid: the energy source for the upper layer.

## Anti-patterns
- **Newtonian substitution for mucus**: a wrong elastic character leads to total dissipation of ciliary energy.
- **Mucus too stiff**: cilia cannot penetrate; transport inhibited.
- **Stokes-approximation overreach**: the non-Newtonian zero-shear-rate Reynolds number is not the only governing parameter (Astarita & Marucci [1]); and rheological non-linearity "would annihilate the advantage of the linearization".
- **Clinical framing**: clearance deficiency (cystic fibrosis, immotile cilia syndrome) causes recurrent bronchitis/pneumonia, chronic otitis media, conductive hearing loss (Corkey et al. [7], Wasserman [24]).

## Key Takeaways
1. One unified two-layer model can serve all ciliated epithelia because structure and beat features are common.
2. Penetration during the effective stroke is experimentally observed and required for propulsion.
3. No prior model class handles both the two-layer structure and mucus non-linearity; that gap defines the thesis.
4. Scope (p. 10): infinite plane double layer, Fourier series in x₂ and t, flat surfaces first, numerical computation, explicit Newtonian comparison, generalized functions throughout.
5. Parts of this work were previously published as Liron & Rozenson [13] and Rozenson & Liron [18] (p. 11).

## Connects To
- **Ch 3**: justifies the force-pulse postulate and averaging announced here.
- **Blake [2]**: the lower-layer solution is inherited from his work; the deep-penetration assumption is rejected.
