from fractions import Fraction
from itertools import permutations
import json
from math import factorial, isqrt
from pathlib import Path


def prime(number):
    return number >= 2 and all(number % divisor for divisor in range(2, isqrt(number) + 1))


def check_siblings(terms, values, factorials, require):
    counts = {'primorial_divisibility': 0, 'prime_tail': 0, 'prime_period': 0,
              'arbitrary_weight_tails': 0, 'quet_constant_tail': 0,
              'companion_coefficients': 0, 'companion_prefix_terms': 0,
              'partition_prime_characterization': 0, 'asymptotic_bounds': 0,
              'hypertree_prime_formula': 0, 'matrix_orbit_shapes': 0,
              'rectangular_band_counts': 0}
    primes = [number for number in range(2, len(terms)) if prime(number)]
    for index in range(1, len(terms)):
        for divisor, partner, quotient in terms[index]:
            if divisor >= 2 and partner >= 2:
                for modulus in primes:
                    if modulus > index:
                        break
                    require(quotient % modulus == 0, 'primorial divisibility')
                    counts['primorial_divisibility'] += 1
    for modulus in primes:
        if modulus > 64:
            break
        for index in range(modulus, 257):
            for weights, sequence in values.items():
                left_weight, right_weight = weights
                require((sequence[index] - left_weight * right_weight ** index
                         - right_weight * left_weight ** index) % modulus == 0,
                        'prime-modulus tail')
                require((sequence[index + modulus - 1] - sequence[index]) % modulus == 0,
                        'prime-modulus period')
                counts['prime_tail'] += 1
                counts['prime_period'] += 1
    for modulus in range(1, 65):
        onset = modulus if prime(modulus) else 2 * modulus + 4
        for index in range(onset, 257):
            for weight in [lambda divisor, partner: divisor - partner,
                           lambda divisor, partner: (-1) ** divisor * (divisor + partner ** 2),
                           lambda divisor, partner: divisor ** 3 - 7 * partner ** 2]:
                total = sum(weight(divisor, partner) * quotient for divisor, partner, quotient in terms[index])
                endpoint = weight(1, index) + weight(index, 1)
                require((total - endpoint) % modulus == 0, 'integer-weight tail')
                counts['arbitrary_weight_tails'] += 1
        for index in range(2 * modulus, 257):
            total = sum(factorials[index] // factorials[divisor] for divisor, partner, quotient in terms[index])
            require((total - 1) % modulus == 0, 'A057625 tail')
            counts['quet_constant_tail'] += 1
    companion = {entry: {} for entry in ['A038041', 'A057625', 'A061095', 'A236696', 'A320444']}
    for index in range(1, 65):
        partition_sum = Fraction(0)
        boxes_sum = Fraction(0)
        forest_sum = Fraction(0)
        for divisor, partner, quotient in terms[index]:
            partition_sum += Fraction(factorials[index], factorials[partner] * factorials[divisor] ** partner)
            boxes_sum += Fraction(factorials[index], factorials[divisor] ** partner)
            forest_sum += factorials[index] * Fraction(divisor ** (divisor - 1), factorials[divisor]) ** partner / factorials[partner]
        for total in [partition_sum, boxes_sum, forest_sum]:
            require(total.denominator == 1, 'integral companion count')
        companion['A038041'][index] = int(partition_sum)
        companion['A061095'][index] = int(boxes_sum)
        companion['A236696'][index] = int(forest_sum)
        companion['A057625'][index] = sum(factorials[index] // factorials[divisor] for divisor, partner, quotient in terms[index])
        require((partition_sum == 2) == prime(index), 'partition prime criterion')
        counts['partition_prime_characterization'] += 1
        error = boxes_sum / factorials[index] - 1
        require(error >= 0 and error ** 2 * 2 ** index <= index ** 2, 'boxes asymptotic bound')
        counts['asymptotic_bounds'] += 1
        hypertree = sum((Fraction(quotient) * Fraction(index + 1, divisor) ** (partner - 1)
                         for divisor, partner, quotient in terms[index]), Fraction(0))
        require(hypertree.denominator == 1, 'integral hypertree count')
        companion['A320444'][index + 1] = int(hypertree)
        if prime(index):
            require(hypertree == 1 + (index + 1) ** (index - 1), 'corrected hypertree prime index')
            counts['hypertree_prime_formula'] += 1
        if index <= 24:
            coefficients = {entry: Fraction(0) for entry in ['A038041', 'A057625', 'A061095', 'A236696']}
            for outer in range(1, index + 1):
                for inner in range(1, index // outer + 1):
                    if outer * inner != index:
                        continue
                    coefficients['A038041'] += Fraction(1, factorial(inner) * factorial(outer) ** inner)
                    coefficients['A057625'] += Fraction(1, factorial(outer))
                    coefficients['A061095'] += Fraction(1, factorial(outer) ** inner)
                    coefficients['A236696'] += Fraction(outer ** (outer - 1), factorial(outer)) ** inner / factorial(inner)
            for entry, coefficient in coefficients.items():
                require(coefficient * factorials[index] == companion[entry][index], 'companion coefficient')
                counts['companion_coefficients'] += 1
    companion['A320444'][0] = companion['A320444'][1] = 1
    prefixes = json.loads(Path(__file__).with_name('sibling-prefixes.json').read_text())
    for entry, record in prefixes.items():
        for index, expected in record['terms']:
            require(companion[entry][index] == int(expected), 'companion OEIS prefix')
            counts['companion_prefix_terms'] += 1
    for index in range(1, 6):
        matrix_count = 0
        operations = set()
        for rows, columns, quotient in terms[index]:
            representatives = set()
            for array in permutations(range(index)):
                orbit = []
                for row_order in permutations(range(rows)):
                    for column_order in permutations(range(columns)):
                        orbit.append(tuple(array[row * columns + column] for row in row_order for column in column_order))
                representatives.add(min(orbit))
                positions = {entry: position for position, entry in enumerate(array)}
                operations.add(tuple(array[(positions[left] // columns) * columns + positions[right] % columns]
                                     for left in range(index) for right in range(index)))
            require(len(representatives) == quotient, 'matrix orbit count')
            matrix_count += len(representatives)
            counts['matrix_orbit_shapes'] += 1
        require(matrix_count == values[1, 1][index] == len(operations), 'rectangular-band count')
        counts['rectangular_band_counts'] += 1
    residues = [values[1, 1][index] % 89 for index in range(41, 45)]
    require(residues == [2, 2, 2, 76], 'three-residue counterexample')
    for index, divisors, expected in [(42, [2, 3, 6], [82, 71, 25]), (44, [2, 4], [21, 16])]:
        require([factorials[index] // (factorials[divisor] * factorials[index // divisor]) % 89
                 for divisor in divisors] == expected, 'counterexample quotient residues')
    large_index = 468
    large_value = sum(factorial(large_index) // (factorial(divisor) * factorial(large_index // divisor))
                      for divisor in range(1, large_index + 1) if large_index % divisor == 0)
    require(len(str(large_value)) == 1007, 'A121860 decimal length at n=468')
    return {'counts': counts, 'prime_moduli': [2, 61], 'prime_base_indices': 'p <= n <= 256',
            'primorial_indices': [4, len(terms) - 1], 'companion_indices': [1, 64],
            'companion_coefficient_indices': [1, 24], 'matrix_and_band_indices': [1, 5],
            'integer_weight_functions': 3, 'counterexample': {'modulus': 89, 'indices': [41, 44], 'residues': residues},
            'corrected_hypertree_example': {'H(3)': 4, 'H(4)': 17},
            'A121860_decimal_length': {'index': large_index, 'digits': len(str(large_value))},
            'scope': 'Exact finite checks of the stated ranges; infinite statements have written proofs.'}
