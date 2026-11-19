# Common permutation laws for rank-one unitary perturbations

This research note studies one probability law on permutations that reproduces the squared moduli of all minors of a fixed unitary matrix. The matrix acts by a phase $\zeta$ on the line spanned by a unit vector $v$ and fixes its orthogonal complement. Put $p_i=|v_i|^2$ and $c=|\zeta-1|^2$.

In dimension five the law exists exactly when

$$c\left(1-\sum_i p_i^2\right)\le3.$$

For a Householder reflection this becomes $\sum_i |v_i|^4\ge1/4$. In arbitrary dimension the note gives a signed-partition criterion, an induction proof that produces an optimal separating metric, and a fractional-triangle-packing construction of a law whenever the criterion holds. In dimension six the example $p=(5,3,1,1,1,1)/12$ requires nonidentity mass $109/108$ for a reflection, although the expected number of moved points is $53/18<3$.

[Read the paper](paper.pdf) · [LaTeX source](paper.tex)

## Reproduce the calculations

Use Python 3.10 or newer with its standard library. Run the commands from this directory.

~~~sh
python verify.py
python construct_law.py --p 2/5 3/20 3/20 3/20 3/20 --c 4
python construct_law.py --p 5/12 3/12 1/12 1/12 1/12 1/12 --c 4
~~~

The constructor uses exact rational simplex arithmetic with Bland's rule. It returns a probability for each nontrivial cycle and the identity, or a rational separating edge vector. Labels in the JSON start at one. The verification script compares the program with direct enumeration of signed partitions, evaluates the historical five-point separator on all 120 permutations, and independently expands determinants over Gaussian integers to check each minor event in three examples.

## Prior work

Drury's Corollary 4 (p. 161) proves the ordinary determinantal inclusion in every dimension for the entire rank-one phase family. Those convex weights can depend on the spectra. His Theorem 7 (pp. 163–164) gives a common minor law for reflections in dimensions at most four. In five dimensions his separator (p. 163) rules out the uniform reflection and gives the necessary bound above when evaluated at general masses. His expansion (7) (p. 162) directly gives a common law whenever $c(1-\sum_i p_i^2)\le2$. [Drury, original paper](https://doi.org/10.1016/S0024-3795%2899%2900173-1)

This note proves sufficiency up to the sharp five-point bound and constructs an optimizer for product edge weights in every dimension. At $p=(8,3,3,3,3)/20$ and $c=4$, twenty oriented three-cycles give an exact law. Even the single-entry probabilities of this example lie outside the convex hull of coordinate embeddings of the four-point phase family. The proof and two independent ways to check the example are in Sections 2 and 5 and in the verification script.

Laurent classified the five-point metric vertices, including the constant $2/3$ vertex and its switchings. The induction in this note proves that a product objective on the triangle-perimeter relaxation admits an optimizer of signed-partition form. Section 6 compares the result with the existing metric, cut and triangle-packing literature. [Laurent, §4.3](https://pure.uvt.nl/ws/portalfiles/portal/1215595/Metric.pdf) · [Bennett, English and Talanda-Fisher, Prop. 2.1](https://arxiv.org/abs/1809.05028) · [del Valle, Dukes and Garaschuk](https://arxiv.org/abs/1608.06017) · [Colbourn, Lemma 2.1](https://combinatorialpress.com/article/jcmcc/Volume%2014/vol-14-paper%2012.pdf)
