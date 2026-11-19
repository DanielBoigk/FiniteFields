import Mathlib.Algebra.Field.ZMod

namespace FiniteFields

/-- Mathlib already provides the field structure on `ZMod p` for prime `p`. -/
example (p : ℕ) [Fact p.Prime] : Field (ZMod p) := inferInstance

end FiniteFields
