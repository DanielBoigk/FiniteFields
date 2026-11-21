/-
The rank encoding of ADR-0005 is order-preserving (backlog R2).

Below `ω^(ω^ω)`, an element of Conway's field `On₂` is a finite set of "digits" (DiMuro,
*On On_p*, arXiv:1108.0962v3, Lemma 1.1), nim-addition is the symmetric difference of digit
sets, and the ordinal order is the colexicographic order: compare the largest digit in which two
elements differ. A subfield `𝔽_{2^n}` is an `𝔽₂`-subspace; take its reduced echelon basis
`b 0, …, b (n-1)`: the largest digits `p i = max (b i)` increase strictly with `i`, and `p i`
occurs in no other `b j`. The element with rank `k = ∑_{i ∈ e} 2^i` is the symmetric difference
`spanSum b e` of the `b i`, `i ∈ e`. This file proves that this map is strictly monotone, from the
binary order of ranks (equivalently colex on index sets) to colex on digit sets.

Only this combinatorial step is formalized. The identification of the ordinal order on `On₂`
with colex on digit sets is DiMuro's Lemma 1.1 (see the vault note *Canonical extension of the
nimbers*), not formalized here.
-/
import Mathlib.Combinatorics.Colex

namespace FiniteFields.RankEncoding

open Finset Finset.Colex
open scoped symmDiff

variable {α ι : Type*} [LinearOrder α] [LinearOrder ι]

/-- The sum over `𝔽₂` (the symmetric difference) of the digit sets `b j`, `j ∈ e`: a digit belongs
to it iff it lies in an odd number of them. -/
def spanSum (b : ι → Finset α) (e : Finset ι) : Finset α :=
  (e.biUnion b).filter fun x => Odd (e.filter fun j => x ∈ b j).card

omit [LinearOrder ι] in
theorem mem_spanSum {b : ι → Finset α} {e : Finset ι} {x : α} :
    x ∈ spanSum b e ↔ Odd (e.filter fun j => x ∈ b j).card := by
  unfold spanSum
  constructor
  · exact fun h => (mem_filter.1 h).2
  · intro h
    refine mem_filter.2 ⟨?_, h⟩
    have hpos : 0 < (e.filter fun j => x ∈ b j).card := by
      obtain ⟨k, hk⟩ := h
      omega
    obtain ⟨j, hj⟩ := card_pos.1 hpos
    exact mem_biUnion.2 ⟨j, (mem_filter.1 hj).1, (mem_filter.1 hj).2⟩

/-- Adding one more basis vector adds its digit set (symmetric difference). -/
theorem spanSum_insert {b : ι → Finset α} {e : Finset ι} {j : ι} (hj : j ∉ e) :
    spanSum b (insert j e) = b j ∆ spanSum b e := by
  ext x
  rw [mem_symmDiff, mem_spanSum, mem_spanSum, filter_insert]
  by_cases hx : x ∈ b j
  · simp only [hx, ↓reduceIte]
    rw [card_insert_of_notMem (by simp [hj])]
    simp [Nat.odd_add_one]
  · simp [hx]

/-- **The rank encoding is order-preserving.** For a reduced echelon family `b` (largest digits
strictly increasing, each largest digit absent from the other members), `e ↦ spanSum b e` is
strictly monotone for the colexicographic orders. -/
theorem spanSum_strictMono (b : ι → Finset α) (hne : ∀ i, (b i).Nonempty)
    (hp : StrictMono fun i => (b i).max' (hne i))
    (hred : ∀ i j, i ≠ j → (b i).max' (hne i) ∉ b j) :
    StrictMono fun e : Colex (Finset ι) => toColex (spanSum b (ofColex e)) := by
  intro e e' h
  obtain ⟨hee, hmax⟩ := Colex.lt_iff_max'_mem.1 h
  set s := ofColex e
  set t := ofColex e'
  have hst : s ≠ t := fun h' => hee (by simpa [s, t] using congrArg toColex h')
  set i := (s ∆ t).max' (symmDiff_nonempty.2 hst)
  have hi_t : i ∈ t := hmax
  have hi_s : i ∉ s := by
    have := max'_mem (s ∆ t) (symmDiff_nonempty.2 hst)
    rw [mem_symmDiff] at this
    tauto
  -- above `i` the index sets agree
  have hagree : ∀ j, i < j → (j ∈ s ↔ j ∈ t) := by
    intro j hij
    by_contra hne'
    have hj : j ∈ s ∆ t := by rw [mem_symmDiff]; tauto
    exact absurd (le_max' _ j hj) (not_le.2 hij)
  set m := (b i).max' (hne i)
  -- a digit above `m` only occurs in `b j` with `j > i`
  have habove : ∀ x j, m < x → x ∈ b j → i < j := by
    intro x j hmx hxj
    exact hp.lt_iff_lt.1 (lt_of_lt_of_le hmx (le_max' _ x hxj))
  have hfilter : ∀ x, m < x → (s.filter fun j => x ∈ b j) = (t.filter fun j => x ∈ b j) := by
    intro x hmx
    ext j
    simp only [mem_filter]
    constructor
    · rintro ⟨hjs, hxj⟩; exact ⟨(hagree j (habove x j hmx hxj)).1 hjs, hxj⟩
    · rintro ⟨hjt, hxj⟩; exact ⟨(hagree j (habove x j hmx hxj)).2 hjt, hxj⟩
  -- `m` lies in `b j` only for `j = i`
  have honly : ∀ (u : Finset ι), (u.filter fun j => m ∈ b j) = u.filter (· = i) := by
    intro u
    ext j
    simp only [mem_filter]
    constructor
    · rintro ⟨hju, hmj⟩
      refine ⟨hju, ?_⟩
      by_contra hji
      exact hred i j (Ne.symm hji) hmj
    · rintro ⟨hju, rfl⟩
      exact ⟨hju, max'_mem _ _⟩
  have hm_t : m ∈ spanSum b t := by
    rw [mem_spanSum, honly, filter_eq' t i]
    simp [hi_t]
  have hm_s : m ∉ spanSum b s := by
    rw [mem_spanSum, honly, filter_eq' s i]
    simp [hi_s]
  show toColex (spanSum b s) < toColex (spanSum b t)
  rw [toColex_lt_toColex_iff_exists_forall_lt]
  refine ⟨m, hm_t, hm_s, fun x hxs hxt => ?_⟩
  rcases lt_trichotomy x m with hlt | heq | hgt
  · exact hlt
  · exact absurd (heq ▸ hxs) hm_s
  · exact absurd (by rw [mem_spanSum, ← hfilter x hgt]; exact mem_spanSum.1 hxs) hxt

/-- The same, with index sets ordered by the rank `∑_{i ∈ e} 2^i` (binary numbers): the element
of larger rank is the larger ordinal. -/
theorem spanSum_lt_of_rank_lt (b : ℕ → Finset α) (hne : ∀ i, (b i).Nonempty)
    (hp : StrictMono fun i => (b i).max' (hne i))
    (hred : ∀ i j, i ≠ j → (b i).max' (hne i) ∉ b j) {e e' : Finset ℕ}
    (h : ∑ i ∈ e, 2 ^ i < ∑ i ∈ e', 2 ^ i) :
    toColex (spanSum b e) < toColex (spanSum b e') :=
  spanSum_strictMono b hne hp hred
    ((geomSum_lt_geomSum_iff_toColex_lt_toColex le_rfl).1 h)

end FiniteFields.RankEncoding
