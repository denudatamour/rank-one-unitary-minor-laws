"""Check the examples by signed partitions, permutation events and exact determinants."""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, prod

from construct_law import construct


def partitions(count):
    """Restricted growth strings enumerate set partitions once."""
    def extend(labels):
        if len(labels) == count:
            yield tuple(labels)
        else:
            for block in range(max(labels) + 2):
                yield from extend(labels + [block])

    yield from extend([0])


def partition_optimum(sizes):
    total = sum(sizes)
    best = None
    count = 0
    for labels in partitions(len(sizes)):
        for tail in product((-1, 1), repeat=len(sizes) - 1):
            signs = (1,) + tail
            imbalances = [0] * (max(labels) + 1)
            for label, sign, size in zip(labels, signs, sizes):
                imbalances[label] += sign * size
            numerator = 3 * total**2 + sum(imbalances)**2 - 4 * sum(x*x for x in imbalances)
            best = numerator if best is None else max(best, numerator)
            count += 1
    return Fraction(best, 12 * total**2), count


def gaussian_determinant(matrix):
    """Leibniz expansion over Gaussian integers, including the empty determinant."""
    order = len(matrix)
    real = imag = 0
    for sigma in permutations(range(order)):
        inversions = sum(sigma[i] > sigma[j] for i, j in combinations(range(order), 2))
        x, y = (-1)**inversions, 0
        for i, j in enumerate(sigma):
            a, b = matrix[i][j]
            x, y = x*a - y*b, x*b + y*a
        real += x
        imag += y
    return real, imag


def direct_minor(sizes, phase, rows, cols):
    """Remove fixed zero coordinates, then use real diagonal similarity to clear roots."""
    u, v, d = phase
    scale = d * sum(sizes)
    if {i for i in rows if sizes[i] == 0} != {j for j in cols if sizes[j] == 0}:
        return Fraction(0)
    rows = [i for i in rows if sizes[i]]
    cols = [j for j in cols if sizes[j]]
    matrix = [[(scale * (i == j) + (u-d) * sizes[j], v * sizes[j])
               for j in cols] for i in rows]
    a, b = gaussian_determinant(matrix)
    return Fraction((a*a+b*b) * prod(sizes[i] for i in rows),
                    scale**(2*len(cols)) * prod(sizes[j] for j in cols))


def explicit_law(data):
    n = len(data["p"])
    atoms = [(tuple(range(n)), Fraction(data["identity"]))]
    for entry in data["cycles"]:
        cycle = [label - 1 for label in entry["cycle"]]
        sigma = list(range(n))
        for i, j in zip(cycle, cycle[1:] + cycle[:1]):
            sigma[i] = j
        atoms.append((tuple(sigma), Fraction(entry["probability"])))
    assert all(probability >= 0 for _, probability in atoms)
    assert sum(probability for _, probability in atoms) == 1
    return atoms


def check_minor_events(sizes, phase):
    u, v, d = phase
    assert u*u + v*v == d*d
    c = Fraction(2*(d-u), d)
    p = [Fraction(x, sum(sizes)) for x in sizes]
    data = construct(p, c)
    assert data["admissible"]
    atoms = explicit_law(data)
    checked = 0
    for k in range(len(p)+1):
        sets = list(combinations(range(len(p)), k))
        for cols in sets:
            distribution = {rows: Fraction(0) for rows in sets}
            for sigma, mass in atoms:
                distribution[tuple(sorted(sigma[j] for j in cols))] += mass
            for rows, probability in distribution.items():
                assert direct_minor(sizes, phase, rows, cols) == probability
                checked += 1
    assert checked == comb(2*len(p), len(p))
    return checked


def check_old_witness():
    """Integer units of 1/6 for Drury's two-set witness (1999, printed p. 163)."""
    pairs = list(map(frozenset, combinations(range(5), 2)))
    observed = set()
    for sigma in permutations(range(5)):
        score = 0
        for subset in pairs:
            image = frozenset(sigma[i] for i in subset)
            distance = len(subset - image)
            score += (3, -1, 1)[distance]
        observed.add(Fraction(score, 6))
    assert observed == {-1, 0, 1, 5}
    return min(observed)


def check_cubic_weights(p, c):
    """Check the signed cubic formula even when its identity mass is negative."""
    triangles = {triple: c * prod(p[i] for i in triple)
                 for triple in combinations(range(len(p)), 3)}
    transpositions = {(i, j): c * p[i] * p[j] * (p[i] + p[j])
                      for i, j in combinations(range(len(p)), 2)}
    identity = 1 - c * (1 - sum(x**3 for x in p)) / 3
    assert identity + sum(transpositions.values()) + 2 * sum(triangles.values()) == 1
    for (i, j), mass in transpositions.items():
        triangle_mass = sum(w for triple, w in triangles.items() if i in triple and j in triple)
        assert mass + triangle_mass == c * p[i] * p[j]
    return identity


def main():
    cases = [
        ((1,), Fraction(0)),
        ((0, 1, 1, 1, 1), Fraction(1, 4)),
        ((1, 1, 1, 1, 1), Fraction(4, 15)),
        ((8, 3, 3, 3, 3), Fraction(1, 4)),
        ((5, 3, 1, 1, 1, 1), Fraction(109, 432)),
        ((6, 4, 3, 2, 1, 1, 1), None),
    ]
    for sizes, expected in cases:
        optimum, count = partition_optimum(sizes)
        p = [Fraction(x, sum(sizes)) for x in sizes]
        data = construct(p, Fraction(4))
        assert optimum == Fraction(data["F"])
        if expected is not None:
            assert optimum == expected
        print(f"sizes {sizes}, {count} signed partitions, F = {optimum}")

    boundary_p = [Fraction(x, 20) for x in (8, 3, 3, 3, 3)]
    h = 4 * (1 - sum(x*x for x in boundary_p))
    assert h == 3 and 1 - h/2 == Fraction(-1, 2)
    diagonal_entry = (1 - 2 * boundary_p[0])**2
    assert diagonal_entry == Fraction(1, 25) < Fraction(1, 4)
    assert sum(x**3 for x in boundary_p) == Fraction(31, 400)
    assert check_cubic_weights(boundary_p, Fraction(4)) == Fraction(-23, 100)
    assert check_cubic_weights(boundary_p, Fraction(3)) == Fraction(31, 400)
    assert check_cubic_weights([Fraction(1, 2), Fraction(1, 2), Fraction(0)], Fraction(4)) == 0
    boundary = construct(boundary_p, Fraction(4))
    assert boundary["identity"] == "0" and len(boundary["cycles"]) == 20
    for atom in boundary["cycles"]:
        assert Fraction(atom["probability"]) == (Fraction(2, 25) if 1 in atom["cycle"]
                                                  else Fraction(1, 200))
    sizes = (5, 3, 1, 1, 1, 1)
    p = [Fraction(x, 12) for x in sizes]
    outside = construct(p, Fraction(4))
    assert 4*(1-sum(x*x for x in p)) == Fraction(53, 18)
    assert not outside["admissible"]
    assert Fraction(outside["necessary_nonidentity_mass"]) == Fraction(109, 108)
    assert 3*Fraction(outside["F"])*12**2 - 2*53 == 3
    # Each weight applies to a group triple in the six-partite graph.
    explicit_packing = {}
    for triple in combinations(range(6), 3):
        groups = (0 in triple, 1 in triple)
        explicit_packing[triple] = {(True, True): 3, (True, False): Fraction(2, 3),
                                    (False, True): 0, (False, False): Fraction(1, 6)}[groups]
    assert sum(explicit_packing.values()) == Fraction(50, 3)
    for i, j in combinations(range(6), 2):
        used = sum(x for triple, x in explicit_packing.items() if i in triple and j in triple)
        assert sizes[i]*sizes[j] - used == (3 if (i, j) == (0, 1) else 0)
    w = {tuple(atom["edge"]): Fraction(atom["value"])
         for atom in outside["separating_edge_weights"]}
    assert all(0 <= weight <= 1 for weight in w.values())
    assert all(sum(w[edge] for edge in combinations(triple, 2)) <= 2
               for triple in combinations(range(1, 7), 3))
    assert 4*sum(p[i-1]*p[j-1]*weight for (i, j), weight in w.items()) == Fraction(109, 108)

    print(f"minimum of the published two-set witness on all 120 permutations = {check_old_witness()}")
    for sizes, phase in [((8, 3, 3, 3, 3), (-1, 0, 1)),
                         ((0, 3, 2, 1, 1), (-3, 4, 5)),
                         ((5, 3, 1, 1, 1, 1), (0, 1, 1))]:
        print(f"sizes {sizes}, phase {phase}: {check_minor_events(sizes, phase)} minor events verified")


if __name__ == "__main__":
    main()
