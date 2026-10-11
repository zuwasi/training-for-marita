# Appendices A-B and References (printed pp. 68-73, PNG 70-75; figure plates PNG 76-81; Hebrew back matter PNG 82-84)

## Core Idea
Appendix A proves the stress tensor has only τ₂₃ = τ₃₂; Appendix B proves the O(ε) wavy-boundary corrections leave all mean values unchanged; the 26-entry reference list anchors the model landscape.

## Appendix A: Form of the Stress Tensor (pp. 68-69)
Starting from the linearized Oldroyd equation (4.4) for the components τ₁₁, τ₁₂ = τ₂₁, τ₂₂, τ₃₃, τ₁₃ = τ₃₁ (all with the same symbol τ, τ = τ(x₂,x₃,t), v₂ = v₂(x₃,t)):
$$\tau+\lambda_1\frac{\partial\tau}{\partial t}+\lambda_1\frac{\partial\tau}{\partial x_2}v_2=0$$
Separate variables τ = T(x₃,t)·X(x₃,x₂):
- $\frac{1}{v_2}+\frac{\lambda_1}{Tv_2}\frac{\partial T}{\partial t}=-\frac{\lambda_1}{X}\frac{\partial X}{\partial x_2}=k(x_3)$
- ODEs: (i) $\partial_tT - T\frac{1}{\lambda}(kv_2-1)=0$; (ii) $\partial_{x_2}X+\frac{1}{\lambda_1}kX=0$
- Solutions as printed: $T=C_T(x_3)\exp[(x_2-x_{20})\frac{k}{\lambda_1}-\frac{t-t_0}{\lambda_1}]$, $X=C_X(x_3)\exp[-k\frac{x_2-x_{20}}{\lambda_1}]$, so $\tau=C(x_3)\exp(\frac{t-t_0}{\lambda_1})$ with C = C_T·C_X.
- **Periodicity forces C(x₃) = 0**: no periodic non-trivial solution exists, so these components vanish. Only τ₂₃ = τ₃₂ survives.

## Appendix B: Vanishing Mean O(ε) Corrections (pp. 70-71)
For the wavy-boundary first-order systems (4.50-4.55) and (4.56-4.61), the (m0) Fourier coefficients of v₃₀, v₂₁, τ₂₃₁, τ₂₂₀, τ₃₃₀ all vanish:
- Key identity: $(v_{20}v_{30})_{m0}=\sum_{m',s}(v_{30})_{m's}(v_{20})_{m+m',s}=0$ because (v₂₀)_{ms} = 0 for s ≠ 0 and (v₃₀)_{m'0} = 0 for s = 0.
- Supporting facts: d(v₃₀)_{m0}/dx₃ = 0; (v₃₀^{(N)})_{m0} = 0 by no-slip at the bottom; (v₃₀^{(NN)})_{m0} = 0 by continuity at the interface; (∂f/∂x₂)v₂₀ terms vanish by periodic integration by parts.
- **Consequence**: the mean x₂-velocity is unchanged at O(ε); the plane penetration solution (4.50) is the full mean answer. The (ms) coefficients are "much more complicated ... not very instructive" and were not solved; the τ₂₃₀ pulse sequence is expected to persist, possibly plus an additive τ₂₃₁ component.
- Procedure for (4.56-4.61): differentiate (4.58), (4.60) in x₃ and (4.57), (4.61) in x₂, then subtract.

## References (26 entries, printed pp. 71-73)
1. Astarita, G. & Marucci, G. (1974), Principles of non-Newtonian Fluid Mechanics.
2. Blake, J.R. (1984), Mechanics of Muco-Ciliary Transport, IMA J. Appl. Math. 32, 69-87 (also Preprint No. 8/83, Australia).
3. Blake, J.R. (1984), Mechanics of Muco-Ciliary Transport, Report to the Fifth Int. Conf. on Physico-Chemical Hydrodynamics, Tel Aviv.
4. Blake, J.R. & Winet, H. (1980), On the Mechanics of Muco-Ciliary Transport, Biorheology 17, 125-134.
5. Brennen, Cr. & Winet, H. (1977), Fluid Mechanics of Propulsion by Cilia and Flagella, Ann. Rev. Fluid Mech. 9, 339-398.
6. Cheung, A.T.W. & Jahn, T.L. (1975), Determination of the Movement Pattern of the Epithelial Cilia of Rabbit Trachea..., in Swimming and Flying in Nature, Plenum, vol. 1, 289-299.
7. Corkey, C.W.B., Levinson, H. & Turner, J.A.B. (1981), The Immobile Cilia Syndrome: A Longitudinal Approach, Am. Rev. Respir. Dis. 124(5), 544-548.
8. Gelman, R.A. & Meyer, F.A. (1979), Mucociliary Transference Rate and Mucus: Dependence on Dynamic Storage and Loss Modulus, Am. Rev. Respir. Dis. 120, 553-557.
9. Gilboa, A. & Silberberg, A. (1976), In situ Rheological Characterization of Epithelial Mucus, Biorheology 13, 59-65.
10. Lighthill, M.J. (1964), Introduction to Fourier Analysis and Generalized Functions, Cambridge UP.
11. Liron, N. & Meyer, F. (1980), Fluid Transport in a Thick Layer on Active Ciliated Surface, Biophys. J. 30, 463-472.
12. Liron, N. & Mochon, S. (1976), The Discrete Cilia Approach to Propulsion of Ciliated Micro-organisms, J. Fluid Mech. 75, 593-607.
13. Liron, N. & Rozenson, M. (1983), Muco-Ciliary Transport, J. Submicrosc. Cytol. 15(1), 317-321.
14. Litt, M. (1981), Physicochemical Determinants of Mucociliary Flow, Chest 80, 846-849.
15. Meyer, F.A. & Silberberg, A. (1978), Structure and Function of Mucus, in The Respiratory Tract Mucus (Ciba Found. Symp.), Elsevier/Excerpta Medica.
16. Meyer, F.A. & Silberberg, A. (1980), The Rheology and Molecular Organization of Epithelial Mucus, Biorheology 17, 163-168.
17. Oldroyd, G.B. (1950), On the Formulation of Rheological Equations of State, Proc. Roy. Soc. A200, 523-541.
18. Rozenson, M. & Liron, N. (1983), Muco-Ciliary Transport, in Biological Structures and Coupled Flows, eds. A. Oplatka & M. Balaban, Academic Press / Balaban ISS.
19. Ross, S.M. & Corrsin, S. (1974), Results of an Analytical Model of Mucociliary Pumping, J. Appl. Phys. 37(3), 333-340.
20. Sanderson, M.S. & Sleigh, M.S. (1981), Ciliary Activity of Cultured Rabbit Tracheal Epithelium: Beat Pattern and Metachrony, J. Cell Sci. 47, 331-347.
21. Sleigh, M.A. (1962), The Biology of Cilia and Flagella, Pergamon Press.
22. Spungin, B. & Silberberg, A. (1984), Stimulation of Mucus Secretion, Ciliary Activity, and Transport in Frog Palate Epithelium, Am. J. Physiol. 247.
23. Sugino, K. & Naito, Y. (1982), Simulated Cross-Bridge Patterns Corresponding to Ciliary Beating in Paramecium, Nature 295, 609-611.
24. Wasserman, S.J. (1984), Ciliary Function and Disease, J. Allergy Clin. Immunol. 73(1), 17-19.
25. Winet, H. & Blake, J.R. (1980), On the Mechanics of the Mucociliary Transport. [II] Observations of the channel model, Biorheology 17, 135-150.
26. Yates, G.T., Wu, T.Y., Johnson, R.E., Cheung, A.T.W. & Frand, C.L. (1980A), A Theoretical and Experimental Study of Tracheal Mucociliary Transport, Biorheology 14, 151-162.

## Figure Plates (PNG 76-81; no printed numbers)
All plots: time axis 0-0.090 s (one beat period), depth 0-7.2 μm (free surface at 7.20 μm), plotted quantity v₂/F₀ (m² sec/kg).

| Fig | Quantity | λ₁, λ₂ | μ₁ | μ₂ | Mean V |
|-----|----------|--------|----|----|--------|
| 6 | NN velocity | 5 s, 0.05 s | 0.01 P | 50 P | 0.6 |
| 7 | NN velocity | 5 s, 0.05 s | 0.01 P | 100 P | 0.6 |
| 8 | NN velocity | 5 s, 0.05 s | 0.02 P | 50 P | 0.3 |
| 9 | NN velocity | 5 s, 0.05 s | 0.03 P | 50 P | 0.2 |
| 10 | Tangent stress τ₂₃ | 5 s, 0.05 s | 0.01 P | 50 P | 0.6 |
| 11 | Tangent stress τ₂₃ | 10 s, 0.01 s | 0.01 P | 50 P | 0.6 |

Shape: velocity = sharp spike near the beat, then a flat, slowly decaying depth-uniform plateau. Stress = large negative walls at the period edges (δ-force artifacts) around a near-constant plateau.

## Key Takeaways
1. The τ₂₃-only stress form is not an assumption of convenience: it is proved by periodicity (App. A).
2. Mean transport is robust to boundary waviness at O(ε) (App. B).
3. The published companions to this thesis are [13] and [18] (Liron/Rozenson).

## Connects To
- **Ch 4a**: uses τ₂₃-only form and the vanishing means.
- **Ch 5**: interprets Figs. 6-11.
