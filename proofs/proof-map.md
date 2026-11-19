# Proof map

The statements below refer to the research note in [paper/main.tex](../paper/main.tex).

| Statement | Assumptions and location | Main dependencies |
| --- | --- | --- |
| Theorem 1: common-law criterion | Unit vector v, unit phase zeta; n >= 1; Sections 1-4 | Rank-one minor identities, triangle-packing linear-programming duality, Lemma 2, signed-partition evaluation |
| Lemma 2: metric optimizer | Arbitrary nonnegative masses; Section 3 | Extreme-point alternatives; zero-edge merging; unit-edge peeling; zero and unit duplications |
| Five- and six-point evaluations | Physical phase range 0 <= c <= 4; Section 4 | Signed partitions, cut bounds and the largest two masses |
| Uniform phase threshold | n >= 2; n=1 handled separately | Signed-block Cauchy-Schwarz bound and extremizers |
| Constructions and obstructions | Exact specified masses and phases; Section 5 | Exchange constraints and rational packing witnesses |
| Multipartite graph interpretation | Nonnegative integer group sizes; Section 5 | Aggregating fractional triangle packings by group triples |

Section 6 distinguishes the claims from Drury's ordinary determinantal theorem and earlier metric and design theory. References and priority qualifications remain in the paper.

The [verification programs](../verification/README.md) check the stated examples and arithmetic certificates. The universal statements depend on the proofs and references in the paper.
