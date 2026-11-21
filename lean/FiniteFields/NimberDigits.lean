/-
The rank encoding of ADR-0005 on Conway's nimbers themselves (backlog R2, step 1).

Using CombinatorialGames' `Nimber` (nim-addition and nim-multiplication on all ordinals) we prove
that an ordinal *is* the set of exponents of its binary Cantor normal form, `digits`:
* nim-addition is the symmetric difference of digit sets (`digits_add`);
* the ordinal order is the colexicographic order of digit sets (`digits_strictMono`).
Together with `FiniteFields.RankEncoding.spanSum_lt_of_rank_lt` this gives the order-preservation
of the rank encoding for any reduced echelon family of nimbers (`nimSum_lt_of_rank_lt`), in
particular for the echelon basis of any finite subfield of `On₂`: no appeal to the tower
description of DiMuro's Lemma 1.1 is needed.
-/
import CombinatorialGames.Nimber.SimplestExtension.Basic
import FiniteFields.RankEncoding

namespace FiniteFields.NimberDigits

open Ordinal Nimber Finset Finset.Colex
open scoped symmDiff

noncomputable section
open Classical in
/-- The base-2 digits of an ordinal: the exponents `e` with `2 ^ e` in its binary Cantor normal
form, `o = 2 ^ e₁ + 2 ^ e₂ + ⋯` with `e₁ > e₂ > ⋯`. -/
def digits (o : Ordinal) : Finset Ordinal :=
  if h : o = 0 then ∅ else insert (log 2 o) (digits (o % 2 ^ log 2 o))
termination_by o
decreasing_by exact (mod_lt _ (opow_ne_zero _ two_ne_zero)).trans_le (opow_log_le_self 2 h)
end

theorem digits_zero : digits 0 = ∅ := by
  rw [digits]; simp

theorem digits_of_ne_zero {o : Ordinal} (h : o ≠ 0) :
    digits o = insert (log 2 o) (digits (o % 2 ^ log 2 o)) := by
  rw [digits]; simp [h]

theorem mod_two_opow_lt (o c : Ordinal) : o % 2 ^ c < 2 ^ c :=
  mod_lt _ (opow_ne_zero _ two_ne_zero)

/-- All digits of `o < 2 ^ c` are below `c`. -/
theorem lt_of_mem_digits {o c : Ordinal} (h : o < 2 ^ c) {e : Ordinal} (he : e ∈ digits o) :
    e < c := by
  induction o using WellFoundedLT.induction generalizing c with
  | ind o IH =>
    obtain rfl | ho := eq_or_ne o 0
    · simp [digits_zero] at he
    rw [digits_of_ne_zero ho, mem_insert] at he
    have hl : log 2 o < c := (lt_opow_iff_log_lt one_lt_two ho).1 h
    rcases he with rfl | he
    · exact hl
    · exact (IH _ ((mod_two_opow_lt o _).trans_le (opow_log_le_self 2 ho))
        (mod_two_opow_lt o _) he).trans hl

/-- `2 ^ m + r` with `r < 2 ^ m` has the digits of `r` and `m`. -/
theorem digits_two_opow_add {m r : Ordinal} (hr : r < 2 ^ m) :
    digits (2 ^ m + r) = insert m (digits r) := by
  have hne : 2 ^ m + r ≠ 0 := fun h => opow_ne_zero m two_ne_zero (eq_zero_of_add_right h)
  have hlog : log 2 (2 ^ m + r) = m := by
    simpa using log_opow_mul_add (b := 2) (u := m) (v := 1) (w := r) one_lt_two one_ne_zero hr
  have hmod : (2 ^ m + r) % 2 ^ m = r := by
    simpa [mod_eq_of_lt hr] using mul_add_mod_self (2 ^ m) 1 r
  rw [digits_of_ne_zero hne, hlog, hmod]

/-- The nimber with the digit set `S`: the nim-sum of the `∗2 ^ e`, `e ∈ S`. -/
noncomputable def ofDigits (S : Finset Ordinal) : Nimber := ∑ e ∈ S, ∗(2 ^ e)

theorem ofDigits_lt {S : Finset Ordinal} {c : Ordinal} (h : ∀ e ∈ S, e < c) :
    ofDigits S < ∗(2 ^ c) :=
  (IsGroup.two_opow c).sum_lt fun e he => by
    rw [of.lt_iff_lt]; exact (opow_lt_opow_iff_right one_lt_two).2 (h e he)

/-- The digits of `ofDigits S` are `S`. -/
theorem digits_ofDigits (S : Finset Ordinal) : digits (ofDigits S).val = S := by
  induction S using Finset.induction_on_max with
  | empty => simp [ofDigits, digits_zero]
  | insert a S hlt IH =>
    have ha : a ∉ S := fun h => lt_irrefl a (hlt a h)
    have hS := ofDigits_lt hlt
    have hsum : ofDigits (insert a S) = ∗(2 ^ a) + ofDigits S := by
      simp [ofDigits, sum_insert ha]
    have hval : (ofDigits (insert a S)).val = 2 ^ a + (ofDigits S).val := by
      rw [hsum, ← (IsGroup.two_opow a).add_eq_of_lt hS]; rfl
    rw [hval, digits_two_opow_add (by simpa using val.lt_iff_lt.2 hS), IH]

/-- Every nimber is the nim-sum of the powers of two in its binary Cantor normal form. -/
theorem ofDigits_digits (o : Ordinal) : ofDigits (digits o) = ∗o := by
  induction o using WellFoundedLT.induction with
  | ind o IH =>
    obtain rfl | ho := eq_or_ne o 0
    · simp [digits_zero, ofDigits]
    have hr := mod_two_opow_lt o (log 2 o)
    have hnot : log 2 o ∉ digits (o % 2 ^ log 2 o) := fun h =>
      lt_irrefl _ (lt_of_mem_digits hr h)
    rw [digits_of_ne_zero ho, ofDigits, sum_insert hnot, ← ofDigits,
      IH _ (hr.trans_le (opow_log_le_self 2 ho)), Nimber.two_opow_log_add ho]

theorem ofDigits_add (S T : Finset Ordinal) : ofDigits S + ofDigits T = ofDigits (S ∆ T) := by
  have hsub : S ∩ T ⊆ S ∪ T := inter_subset_union
  have h1 := sum_union_inter (s₁ := S) (s₂ := T) (f := fun e => ∗(2 ^ e))
  have h2 := sum_sdiff hsub (f := fun e => ∗(2 ^ e))
  have hsd : S ∆ T = (S ∪ T) \ (S ∩ T) := symmDiff_eq_sup_sdiff_inf S T
  unfold ofDigits
  rw [hsd, ← h1, ← h2, add_assoc, add_self, add_zero]

theorem ofDigits_digits_val (x : Nimber) : ofDigits (digits x.val) = x := by
  rw [ofDigits_digits]; rfl

/-- **Nim-addition is the symmetric difference of binary digits.** -/
theorem digits_add (x y : Nimber) : digits (x + y).val = digits x.val ∆ digits y.val := by
  calc digits (x + y).val
      = digits (ofDigits (digits x.val) + ofDigits (digits y.val)).val := by
        rw [ofDigits_digits_val, ofDigits_digits_val]
    _ = digits (ofDigits (digits x.val ∆ digits y.val)).val := by rw [ofDigits_add]
    _ = digits x.val ∆ digits y.val := digits_ofDigits _

/-- **The ordinal order is the colexicographic order of binary digits.** -/
theorem digits_strictMono : StrictMono fun o : Ordinal => toColex (digits o) := by
  intro o o' h
  induction o' using WellFoundedLT.induction generalizing o with
  | ind o' IH =>
    have ho' : o' ≠ 0 := (lt_of_le_of_lt zero_le h).ne'
    obtain rfl | ho := eq_or_ne o 0
    · show toColex (digits 0) < toColex (digits o')
      rw [digits_zero, digits_of_ne_zero ho']
      exact toColex_lt_toColex_of_ssubset (empty_ssubset.2 (insert_nonempty _ _))
    have hle : log 2 o ≤ log 2 o' := log_mono_right 2 h.le
    -- every digit of `o` is at most `log 2 o`
    have hdig : ∀ d ∈ digits o, d ≤ log 2 o := fun d hd =>
      Order.lt_succ_iff.1 (lt_of_mem_digits (lt_opow_succ_log_self one_lt_two o) hd)
    rcases hle.lt_or_eq with hlt | heq
    · show toColex (digits o) < toColex (digits o')
      rw [toColex_lt_toColex_iff_exists_forall_lt]
      refine ⟨log 2 o', ?_, fun hm => absurd (hdig _ hm) (not_le.2 hlt),
        fun d hd _ => (hdig d hd).trans_lt hlt⟩
      rw [digits_of_ne_zero ho']
      exact mem_insert_self _ _
    · set l := log 2 o with hl
      have hr := mod_two_opow_lt o l
      have hr' := mod_two_opow_lt o' l
      have hlo : o = 2 ^ l + o % 2 ^ l := (Ordinal.two_opow_log_add ho).symm
      have hlo' : o' = 2 ^ l + o' % 2 ^ l := by rw [heq]; exact (Ordinal.two_opow_log_add ho').symm
      have hrr : o % 2 ^ l < o' % 2 ^ l := by
        have h2 := h
        rw [hlo, hlo'] at h2
        exact (add_lt_add_iff_left _).1 h2
      have IH' := IH _ (hr'.trans_le (by rw [heq]; exact opow_log_le_self 2 ho')) hrr
      have hd : digits o = insert l (digits (o % 2 ^ l)) := digits_of_ne_zero ho
      have hd' : digits o' = insert l (digits (o' % 2 ^ l)) := by
        rw [digits_of_ne_zero ho', ← heq]
      have hcancel : ∀ D : Finset Ordinal, (∀ d ∈ D, d < l) → insert l D \ {l} = D := by
        intro D hD
        ext x
        by_cases hx : x = l
        · subst hx; simpa using fun hm => lt_irrefl _ (hD _ hm)
        · simp [hx]
      show toColex (digits o) < toColex (digits o')
      rw [hd, hd', ← toColex_sdiff_lt_toColex_sdiff (u := {l}) (by simp) (by simp),
        hcancel _ (fun d hd => lt_of_mem_digits hr hd),
        hcancel _ (fun d hd => lt_of_mem_digits hr' hd)]
      exact IH'

/-- The digits of a nim-sum of nimbers are the `𝔽₂`-sum of their digit sets. -/
theorem digits_sum (b : ℕ → Nimber) (e : Finset ℕ) :
    digits (∑ i ∈ e, b i).val = RankEncoding.spanSum (fun i => digits (b i).val) e := by
  induction e using Finset.induction with
  | empty => simp [digits_zero, RankEncoding.spanSum]
  | insert j e hj IH => rw [sum_insert hj, digits_add, IH, RankEncoding.spanSum_insert hj]

/-- **The rank encoding of `On₂` is order-preserving.** Let `b 0, b 1, …` be nimbers forming a
reduced echelon family for the binary digits (the largest digits increase strictly, and the
largest digit of `b i` is not a digit of any other `b j`), e.g. the echelon basis of a finite
subfield of `On₂`. Then the nim-sum `∑_{i ∈ e} b i` increases strictly with the rank
`∑_{i ∈ e} 2 ^ i`: the element of rank `k` is the `k`-th smallest. -/
theorem nimSum_lt_of_rank_lt (b : ℕ → Nimber) (hne : ∀ i, (digits (b i).val).Nonempty)
    (hp : StrictMono fun i => (digits (b i).val).max' (hne i))
    (hred : ∀ i j, i ≠ j → (digits (b i).val).max' (hne i) ∉ digits (b j).val)
    {e e' : Finset ℕ} (h : ∑ i ∈ e, (2 : ℕ) ^ i < ∑ i ∈ e', (2 : ℕ) ^ i) :
    ∑ i ∈ e, b i < ∑ i ∈ e', b i := by
  have key := RankEncoding.spanSum_lt_of_rank_lt (fun i => digits (b i).val) hne hp hred h
  rw [← digits_sum, ← digits_sum] at key
  exact val.lt_iff_lt.1 (digits_strictMono.lt_iff_lt.1 key)

end FiniteFields.NimberDigits
