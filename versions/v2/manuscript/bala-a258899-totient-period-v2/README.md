# Bala's totient-period conjectures: A258899 and A258903

The note proves the conjectures in the COMMENTS of OEIS A258899 (3 August
2025) and A258903 (8 August 2025), and the constant-tail conjecture in
A121860 (8 August 2025). For all integer weights, positive moduli `m`, and
`n >= 2*m+4`, the weighted divisor-factorial sum is congruent to
`u*v**n + v*u**n` modulo `m`. Euler's totient and the Carmichael exponent
are eventual periods. The onset is sufficient, not necessarily least.
For prime moduli `p` the identity and period `p-1` hold from `n=p`.
The added corollaries prove the linked A057625 constant tail and A386877
interior-entry congruences, establish the elementary companion identities,
correct two displayed source formulas, and refute the A121860 three-residue
onset heuristic with residues `2,2,2,76` at indices `41..44` modulo `89`.

## Contents

- `note.tex`, `siblings.tex`: self-contained statement, proof, companion formulas, references,
  and declaration-by-declaration coverage.
- `note.pdf`: typeset note.
- `metadata.json`: bibliographic metadata and covered problem identity.
- `verify.py`, `sibling_checks.py`, `sequence-prefixes.json`, `sibling-prefixes.json`, `checks.json`: exact finite checks,
  source sequence values, and expected output.
- `coverage.json`: statement/source/verification mapping.
- `MANIFEST.sha256`: member hashes for the source package.

## Build

Use Tectonic **0.17.0**:

```sh
tectonic --untrusted note.tex
```

Alternatively, use a standard TeX Live installation:

```sh
pdflatex -halt-on-error note.tex
pdflatex -halt-on-error note.tex
```

The bibliography is included in `note.tex`; no BibTeX database is needed.
The source ZIP contains the TeX source, metadata, checker, sequence prefixes,
coverage table, expected checker output, this guide, and member hashes.
The PDF is supplied separately.

## Exact finite checks

Python 3.11 or newer, standard library only:

```sh
python3 verify.py > recomputed-checks.json
cmp checks.json recomputed-checks.json
sha256sum -c MANIFEST.sha256
```

The checker uses exact integer quotients and independent rational coefficient
extraction. The output lists every finite range and the number of tested
instances. Infinite statements follow from the note's proofs.

For the natural-weight sequences A258899 and A258903, the Lean declaration
`D5.S1.Recurrence.Periodic.BalaDivisorFactorialTotientPeriod.result` states
the two totient-period conclusions with the explicit onset `2*m+4`:

https://github.com/the-omega-institute/trureturing/blob/eb2703740e8b4b97e408a80b2eb0d0940939dfc3/D5/S1/Recurrence/Periodic/BalaDivisorFactorialTotientPeriod.lean

The power lemma `pow_add_totient` is in:

https://github.com/the-omega-institute/trureturing/blob/eb2703740e8b4b97e408a80b2eb0d0940939dfc3/D5/S1/Recurrence/Periodic/StirlingTransformTotientPeriod.lean

The signed-weight extension, factorial strengthening, Carmichael period,
A121860 specialization, and generating-function identities are written proofs.
The prime-modulus onsets, integer-valued weights and additional corollaries
are written proofs. The checker supplies finite controls for every new assertion,
including independent small matrix-orbit and operation-table enumerations.
Devillet and Kiss, arXiv:1806.02073v1 (2018), Proposition 4.13 and
Corollary 4.15, already prove the rectangular-band count and primality criterion.

## Authorship and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang. CC BY 4.0 (`cc-by-4.0`):
https://creativecommons.org/licenses/by/4.0/

Sequence values are attributed to the individual OEIS entries in
`sequence-prefixes.json`; OEIS material is subject to the OEIS End-User License
Agreement, https://oeis.org/wiki/The_OEIS_End-User_License_Agreement.
