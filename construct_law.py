"""Construct a rank-one phase unitary's common minor law over exact rationals.

Example: python construct_law.py --p 2/5 1/5 1/5 1/10 1/10 --c 4
The output lists permutations by their nontrivial cycle, with one-based indices.
A negative answer includes the rational separating edge weights.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path


def triangle_packing(masses, edges, triangles):
    """Primal simplex with Bland's rule, starting with edge slacks as the basis."""
    m, t = len(edges), len(triangles)
    caps = [masses[i] * masses[j] for i, j in edges]
    rows = []
    for e, (i, j) in enumerate(edges):
        row = [F(i in tri and j in tri) for tri in triangles]
        rows.append(row + [F(e == k) for k in range(m)] + [caps[e]])
    basis = list(range(t, t + m))
    gains = [F(1)] * t + [F(0)] * m
    while True:
        entering = next((j for j, gain in enumerate(gains) if gain > 0), None)
        if entering is None:
            break
        ratios = [(row[-1] / row[entering], basis[i], i)
                  for i, row in enumerate(rows) if row[entering] > 0]
        _, _, leaving = min(ratios)
        pivot = rows[leaving][entering]
        pivot_row = [x / pivot for x in rows[leaving]]
        rows[leaving] = pivot_row
        for i, row in enumerate(rows):
            if i != leaving and row[entering]:
                factor = row[entering]
                rows[i] = [x - factor * y for x, y in zip(row, pivot_row)]
        gain = gains[entering]
        gains = [x - gain * y for x, y in zip(gains, pivot_row)]
        basis[leaving] = entering
    solution = [F(0)] * (t + m)
    for i, j in enumerate(basis):
        solution[j] = rows[i][-1]
    packing, slacks = solution[:t], solution[t:]
    cover = [min(F(1), -gains[t + e]) for e in range(m)]
    edge_index = {edge: i for i, edge in enumerate(edges)}
    assert all(x >= 0 for x in solution + cover)
    for tri in triangles:
        assert sum(cover[edge_index[e]] for e in combinations(tri, 2)) >= 1
    for h, (i, j) in enumerate(edges):
        assert slacks[h] + sum(x for tri, x in zip(triangles, packing) if i in tri and j in tri) == caps[h]
    assert sum(packing) == sum(x * y for x, y in zip(caps, cover))
    return packing, slacks, [1 - y for y in cover]


def construct(masses, phase_cost):
    if not masses or any(p < 0 for p in masses) or sum(masses) != 1:
        raise ValueError('The input masses must be nonnegative and sum to one.')
    if not 0 <= phase_cost <= 4:
        raise ValueError('The phase cost must be in [0,4].')
    edges = list(combinations(range(len(masses)), 2))
    triangles = list(combinations(range(len(masses)), 3))
    packing, slacks, dual = triangle_packing(masses, edges, triangles)
    f_value = sum(masses[i] * masses[j] * w for (i, j), w in zip(edges, dual))
    min_mass = phase_cost * f_value
    result = {'p': [str(p) for p in masses], 'phase_cost': str(phase_cost),
              'F': str(f_value), 'necessary_nonidentity_mass': str(min_mass),
              'admissible': min_mass <= 1}
    if min_mass > 1:
        result['separating_edge_weights'] = [
            {'edge': [i + 1, j + 1], 'value': str(w)} for (i, j), w in zip(edges, dual)
        ]
        return result
    cycle_probabilities = []
    for (i, j), slack in zip(edges, slacks):
        if phase_cost * slack:
            cycle_probabilities.append(([i + 1, j + 1], phase_cost * slack))
    for (i, j, k), value in zip(triangles, packing):
        mass = phase_cost * value
        if mass:
            cycle_probabilities.extend([([i + 1, j + 1, k + 1], mass),
                                        ([i + 1, k + 1, j + 1], mass)])
    identity_mass = 1 - min_mass
    assert sum(value for _, value in cycle_probabilities) + identity_mass == 1
    result['identity'] = str(identity_mass)
    result['cycles'] = [{'cycle': cycle, 'probability': str(weight)} for cycle, weight in cycle_probabilities]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p', nargs='+', type=F, required=True, help='squared moduli of the unit vector')
    parser.add_argument('--c', type=F, required=True, help='|zeta-1| squared')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = json.dumps(construct(args.p, args.c), indent=2) + '\n'
    if args.output:
        args.output.write_text(result, encoding='utf-8')
    else:
        print(result, end='')


if __name__ == '__main__':
    main()
