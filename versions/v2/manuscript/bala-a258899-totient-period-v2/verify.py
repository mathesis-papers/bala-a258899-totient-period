import json
from fractions import Fraction
from math import factorial, gcd, lcm
from pathlib import Path
from sibling_checks import check_siblings


def factorization(number):
    factors = []
    prime = 2
    while prime * prime <= number:
        exponent = 0
        while number % prime == 0:
            exponent += 1
            number //= prime
        if exponent:
            factors.append((prime, exponent))
        prime += 1
    if number > 1:
        factors.append((number, 1))
    return factors


def totient(modulus):
    return sum(gcd(residue, modulus) == 1 for residue in range(1, modulus + 1))


def carmichael(modulus):
    result = 1
    for prime, exponent in factorization(modulus):
        component = (prime - 1) * prime ** (exponent - 1)
        if prime == 2 and exponent >= 3:
            component //= 2
        result = lcm(result, component)
    return result


def require(condition, label):
    if not condition:
        raise ValueError(label)


def main():
    maximum_index = 384
    factorials = [factorial(index) for index in range(maximum_index + 1)]
    terms = [[] for index in range(maximum_index + 1)]
    counts = {}
    quotient_count = 0
    for index in range(1, maximum_index + 1):
        for divisor in range(1, index + 1):
            if index % divisor:
                continue
            partner = index // divisor
            quotient, remainder = divmod(factorials[index], factorials[divisor] * factorials[partner])
            require(remainder == 0, 'integral quotient')
            terms[index].append((divisor, partner, quotient))
            quotient_count += 1
    counts['integer_quotients'] = quotient_count
    values = {}
    weights = range(-2, 4)
    for left_weight in weights:
        for right_weight in weights:
            values[left_weight, right_weight] = [1] + [
                sum(left_weight ** divisor * right_weight ** partner * quotient
                    for divisor, partner, quotient in terms[index])
                for index in range(1, maximum_index + 1)
            ]
    coefficient_count = 0
    for left_weight in weights:
        for right_weight in weights:
            for index in range(1, 25):
                coefficient = Fraction(0)
                for outer_index in range(1, index + 1):
                    for inner_index in range(1, index // outer_index + 1):
                        if outer_index * inner_index == index:
                            coefficient += Fraction(left_weight ** outer_index * right_weight ** inner_index,
                                                    factorial(outer_index) * factorial(inner_index))
                require(coefficient * factorial(index) == values[left_weight, right_weight][index],
                        'rational generating-function coefficient')
                coefficient_count += 1
    counts['rational_coefficients'] = coefficient_count
    counts['weight_symmetry'] = 0
    for left_weight in weights:
        for right_weight in weights:
            require(values[left_weight, right_weight] == values[right_weight, left_weight], 'weight symmetry')
            counts['weight_symmetry'] += maximum_index + 1
    for label in ['factorial_tail_terms', 'tail_congruences', 'totient_periods',
                  'carmichael_periods', 'power_congruences', 'constant_tail']:
        counts[label] = 0
    for modulus in range(1, 65):
        phi = totient(modulus)
        exponent = carmichael(modulus)
        require(phi % exponent == 0, 'Carmichael divides totient')
        for base in weights:
            for index in range(modulus, 257):
                for period in [phi, exponent]:
                    require((pow(base, index + period, modulus) - pow(base, index, modulus)) % modulus == 0,
                            'power period')
                    counts['power_congruences'] += 1
        for index in range(2 * modulus + 4, 257):
            for divisor, partner, quotient in terms[index]:
                if divisor >= 2 and partner >= 2:
                    require(quotient % factorials[modulus] == 0, 'factorial strengthening')
                    counts['factorial_tail_terms'] += 1
            require((values[1, 1][index] - 2) % modulus == 0, 'constant tail')
            counts['constant_tail'] += 1
            for left_weight in weights:
                for right_weight in weights:
                    sequence = values[left_weight, right_weight]
                    require((sequence[index] - left_weight * right_weight ** index
                             - right_weight * left_weight ** index) % modulus == 0, 'tail congruence')
                    counts['tail_congruences'] += 1
                    require((sequence[index + phi] - sequence[index]) % modulus == 0, 'totient period')
                    counts['totient_periods'] += 1
                    require((sequence[index + exponent] - sequence[index]) % modulus == 0, 'Carmichael period')
                    counts['carmichael_periods'] += 1
    source_values = json.loads(Path(__file__).with_name('sequence-prefixes.json').read_text())
    counts['oeis_terms'] = 0
    for entry, weights_pair in [('A258899', (2, 1)), ('A258903', (2, 3)), ('A121860', (1, 1))]:
        for index, expected in source_values[entry]['terms']:
            require(values[weights_pair][index] == int(expected), entry + ' prefix')
            counts['oeis_terms'] += 1
    examples = []
    for entry, weights_pair, modulus, onset, period in [
        ('A258899', (2, 1), 11, 11, 10), ('A258903', (2, 3), 7, 7, 6),
    ]:
        sequence = values[weights_pair]
        for index in range(onset, 257):
            require((sequence[index + period] - sequence[index]) % modulus == 0, 'OEIS example')
        examples.append({'entry': entry, 'modulus': modulus, 'onset': onset, 'period': period,
                         'checked_base_indices': [onset, 256]})
    counts['prime_characterization'] = 0
    for index in range(1, 257):
        is_prime = index >= 2 and all(index % divisor for divisor in range(2, int(index ** 0.5) + 1))
        require((values[1, 1][index] == 2) == is_prime, 'prime characterization')
        counts['prime_characterization'] += 1
    sibling_results = check_siblings(terms, values, factorials, require)
    result = {
        'status': 'pass',
        'arithmetic': 'exact integers and fractions',
        'ranges': {'weights': [-2, 3], 'moduli': [1, 64], 'computed_indices': [0, maximum_index],
                   'tail_and_period_base_indices': '2*m+4 <= n <= 256',
                   'power_base_indices': 'm <= n <= 256', 'rational_coefficient_indices': [1, 24]},
        'counts': counts,
        'oeis_examples': examples,
        'sibling_checks': sibling_results,
        'scope': 'Finite checks; the note supplies proofs for all stated parameters.',
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
