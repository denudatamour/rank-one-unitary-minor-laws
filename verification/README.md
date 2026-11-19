# Verification

Use Python 3.10 or newer with the standard library. Run from the repository root, with assertions enabled:

```sh
python -B -S verification/check_hashes.py
python -B -S verify.py
python -B -S construct_law.py --p 2/5 3/20 3/20 3/20 3/20 --c 4
python -B -S construct_law.py --p 5/12 3/12 1/12 1/12 1/12 1/12 --c 4
```

The constructor uses exact rational simplex arithmetic with Bland's rule. It returns cycle probabilities or a rational separating vector; labels in the JSON output start at one.

The verifier independently enumerates signed partitions, evaluates Drury's separator on all 120 five-label permutations, and expands determinants over Gaussian integers for three examples with 252, 252 and 924 minor events. These are finite certificate checks; the universal results depend on the arguments listed in the [proof map](../proofs/proof-map.md).

[Recorded output](exact-checks.txt) · [Environment](environment.json) · [File hashes](SHA256SUMS)

GitHub Actions checks the complete file manifest and reruns the exact verifier. SHA256SUMS covers every tracked file except itself.

## Building the paper

From paper/, run twice to resolve references:

```sh
pdflatex -no-shell-escape -halt-on-error -interaction=nonstopmode -jobname=paper main.tex
```

The root paper.pdf is an identical compatibility copy of paper/paper.pdf; root paper.tex forwards to paper/main.tex. After revising the manuscript, rebuild and inspect the PDF, update the compatibility copy, and regenerate the file manifest.
