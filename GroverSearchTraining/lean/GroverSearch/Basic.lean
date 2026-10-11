import Mathlib

/-!
# Formal companion to Grover, "Quantum Mechanics helps in searching for a needle in a haystack"

L. K. Grover, Phys. Rev. Lett. 79, 325 (1997), arXiv:quant-ph/9706033.

The statements below are written from the normalized claim ledger in
`docs/formalization_map.md`, independently of the Wolfram implementation.
Amplitudes are real throughout, which is exact for Grover's algorithm because every
operator it uses (M, W, R, D, and the phase-π oracle) is a real matrix.
-/

open Matrix Finset

namespace GroverSearch

/-! ## Section 1.1: the single-bit transformation M -/

/-- Grover's `M`, scaled by `√2` so that its entries are `±1` (Sec. 1.1). -/
def Mu : Matrix (Fin 2) (Fin 2) ℝ := !![1, 1; 1, -1]

/-- `M = Mu / √2`. -/
noncomputable def M : Matrix (Fin 2) (Fin 2) ℝ := (Real.sqrt 2)⁻¹ • Mu

/-- Sec. 1.1 / Sec. 5: `M M = I`, so `M` is unitary (it is real and symmetric). -/
theorem M_mul_M : M * M = 1 := by
  have h2 : (Real.sqrt 2)⁻¹ * (Real.sqrt 2)⁻¹ = (1 / 2 : ℝ) := by
    rw [← mul_inv, Real.mul_self_sqrt (by norm_num : (0:ℝ) ≤ 2)]; norm_num
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [M, Mu, Matrix.mul_apply, Fin.sum_univ_two, smul_eq_mul] <;>
    nlinarith [h2]

theorem M_transpose : Mᵀ = M := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [M, Mu]

/-! ## Sections 1.1 and 5: the Walsh-Hadamard transform on `n` bits

Basis states are `n`-bit strings `x : Fin n → Fin 2`. -/

/-- Bitwise dot product `x · y` of two `n`-bit strings. -/
def bitDot {n : ℕ} (x y : Fin n → Fin 2) : ℕ := ∑ k, (x k).val * (y k).val

/-- Unnormalized Walsh-Hadamard matrix `Hu_xy = (-1)^(x·y)`; Grover's `W = 2^(-n/2) Hu`. -/
def Hu (n : ℕ) : Matrix (Fin n → Fin 2) (Fin n → Fin 2) ℝ :=
  fun x y => (-1 : ℝ) ^ bitDot x y

/-- Normalized Walsh-Hadamard matrix `W_ij = 2^(-n/2) (-1)^(i·j)` (Sec. 5). -/
noncomputable def W (n : ℕ) : Matrix (Fin n → Fin 2) (Fin n → Fin 2) ℝ :=
  (Real.sqrt (2 ^ n))⁻¹ • Hu n

lemma inner_bit_sum (a c : Fin 2) :
    ∑ b : Fin 2, (-1 : ℝ) ^ (a.val * b.val) * (-1 : ℝ) ^ (b.val * c.val) =
      if a = c then 2 else 0 := by
  fin_cases a <;> fin_cases c <;> simp [Fin.sum_univ_two] <;> norm_num

/-- Core identity behind `W W = I`: `Hu Hu = 2^n I` for every `n`. -/
theorem Hu_mul_Hu (n : ℕ) : Hu n * Hu n = (2 ^ n : ℝ) • (1 : Matrix _ _ ℝ) := by
  ext x z
  simp only [Matrix.mul_apply, Hu, bitDot, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
  have key : ∀ y : Fin n → Fin 2,
      (-1 : ℝ) ^ (∑ k, (x k).val * (y k).val) * (-1 : ℝ) ^ (∑ k, (y k).val * (z k).val) =
        ∏ k, ((-1 : ℝ) ^ ((x k).val * (y k).val) * (-1 : ℝ) ^ ((y k).val * (z k).val)) := by
    intro y
    rw [Finset.prod_mul_distrib, Finset.prod_pow_eq_pow_sum, Finset.prod_pow_eq_pow_sum]
  simp_rw [key]
  rw [← Fintype.prod_sum (fun k b => (-1 : ℝ) ^ ((x k).val * b.val) * (-1 : ℝ) ^ (b.val * (z k).val))]
  simp_rw [inner_bit_sum]
  by_cases hxz : x = z
  · subst hxz
    simp
  · obtain ⟨k, hk⟩ : ∃ k, x k ≠ z k := by
      by_contra h
      push_neg at h
      exact hxz (funext h)
    rw [if_neg hxz, mul_zero]
    exact Finset.prod_eq_zero (Finset.mem_univ k) (if_neg hk)

lemma sqrt_two_pow_inv_sq (n : ℕ) :
    (Real.sqrt (2 ^ n))⁻¹ * (Real.sqrt (2 ^ n))⁻¹ = ((2 : ℝ) ^ n)⁻¹ := by
  rw [← mul_inv, Real.mul_self_sqrt (by positivity)]

/-- Sec. 5: `W W = I` for every number of bits `n`. -/
theorem W_mul_W (n : ℕ) : W n * W n = 1 := by
  rw [W, Matrix.smul_mul, Matrix.mul_smul, Hu_mul_Hu, smul_smul, smul_smul,
    sqrt_two_pow_inv_sq]
  rw [inv_mul_cancel₀ (by positivity), one_smul]

theorem Hu_transpose (n : ℕ) : (Hu n)ᵀ = Hu n := by
  ext x y
  simp [Hu, bitDot, mul_comm]

theorem W_transpose (n : ℕ) : (W n)ᵀ = W n := by
  rw [W, Matrix.transpose_smul, Hu_transpose]

/-! ## Section 1.1: selective phase rotation preserves probabilities -/

/-- Last paragraph of Sec. 1.1: multiplying an amplitude by `e^{iφ}` leaves
`|amplitude|²` unchanged. -/
theorem phase_rotation_preserves_probability (φ : ℝ) (z : ℂ) :
    Complex.normSq (Complex.exp (φ * Complex.I) * z) = Complex.normSq z := by
  rw [Complex.normSq_mul, Complex.normSq_eq_norm_sq (Complex.exp _),
    Complex.norm_exp_ofReal_mul_I]
  simp

end GroverSearch
