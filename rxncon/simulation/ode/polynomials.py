from typing import List, Callable, Set
from collections import defaultdict
from math import isclose

import typecheck as tc


TRIVIAL_MUL_SYMBOL = '1'
TRIVIAL_ADD_SYMBOL = '0'


### SYMBOLS, i.e. x, y, z ###
class Symbol:
    @tc.typecheck
    def __init__(self, name: str):
        self.name = name

    @tc.typecheck
    def __eq__(self, other: 'Symbol') -> bool:
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def __str__(self) -> str:
        return '{0}'.format(self.name)

    @property
    def is_trivial_mul(self):
        return self.name == TRIVIAL_MUL_SYMBOL

    @property
    def is_trivial_add(self):
        return self.name == TRIVIAL_ADD_SYMBOL


def TrivialMulSymbol():
    return Symbol(TRIVIAL_MUL_SYMBOL)


def TrivialAddSymbol():
    return Symbol(TRIVIAL_ADD_SYMBOL)


### MONOMIALS, i.e. (x^3)*(y^4) ###
class Monomial:
    @tc.typecheck
    def __init__(self, factors: Set['MonomialFactor']):
        self.factors = factors

    @tc.typecheck
    def __eq__(self, other: 'Monomial') -> bool:
        return self.factors == other.factors

    def __hash__(self) -> int:
        return hash('*monomial*{0}'.join(str(x) for x in self.factors))

    @tc.typecheck
    def __mul__(self, other: 'Monomial') -> 'Monomial':
        sym_to_power = defaultdict(int)

        for factor in self.factors:
            sym_to_power[factor.symbol] = factor.power

        for factor in other.factors:
            sym_to_power[factor.symbol] += factor.power

        new_factors = []

        for sym, power in sym_to_power.items():
            if power != 0:
                new_factors.append(MonomialFactor(sym, power))

        if not new_factors:
            return Monomial({MonomialFactor(TrivialMulSymbol(), 0)})

        return Monomial(set(new_factors))

    def __str__(self) -> str:
        return ''.join(str(fac) for fac in self.factors)

class MonomialFactor:
    @tc.typecheck
    def __init__(self, symbol: Symbol, power: int):
        self.symbol = symbol
        self.power = power
        self._validate()

    @tc.typecheck
    def __eq__(self, other: 'MonomialFactor'):
        return self.symbol == other.symbol and self.power == other.power

    def __hash__(self) -> int:
        return hash('*mono-fac*{0}^{1}'.format(str(self.symbol), str(self.power)))

    def __str__(self) -> str:
        if self.power < 0:
            return '{0}^({1})'.format(self.symbol, self.power)
        elif self.power == 0:
            return '1'
        elif self.power > 0:
            return '{0}^{1}'.format(self.symbol, self.power)

    def _validate(self):
        if self.power == 0 and not self.symbol.is_trivial_mul:
            raise ValueError('Only TrivialSymbol may be raised to the zero-th power.')

        elif self.symbol.is_trivial_mul and not self.power == 0:
            raise ValueError('TrivialSymbol may only be raised to the zero-th power.')


def TrivialMonomial():
    return Monomial({MonomialFactor(TrivialMulSymbol(), 0)})


### POLYNOMIALS, i.e. 3*(x^2)*(y^4) - 2*(z^8) ###
class Polynomial:
    @tc.typecheck
    def __init__(self, terms: Set['PolynomialTerm']):
        self.terms = terms
        self._validate()

    @tc.typecheck
    def __add__(self, other: 'Polynomial') -> 'Polynomial':
        new_monomial_to_factor = defaultdict(lambda: 0.0)

        for term in self.terms:
            new_monomial_to_factor[term.monomial] = term.factor

        for term in other.terms:
            new_monomial_to_factor[term.monomial] += term.factor

        new_terms = []

        for new_mon, new_fac in new_monomial_to_factor.items():
            if not isclose(new_fac, 0.0):
                new_terms.append(PolynomialTerm(new_mon, new_fac))

        return Polynomial(set(new_terms))

    @tc.typecheck
    def __eq__(self, other: 'Polynomial') -> bool:
        return self.terms == other.terms

    @tc.typecheck
    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        this_monomial_to_factor = {term.monomial: term.factor for term in self.terms}
        that_monomial_to_factor = {term.monomial: term.factor for term in other.terms}

        new_monomial_to_factor = defaultdict(lambda: 0.0)

        for this_mon, this_fac in this_monomial_to_factor.items():
            for that_mon, that_fac in that_monomial_to_factor.items():
                new_monomial_to_factor[this_mon * that_mon] += this_fac * that_fac

        new_terms = []

        for new_mon, new_fac in new_monomial_to_factor.items():
            if not isclose(new_fac, 0.0):
                new_terms.append(PolynomialTerm(new_mon, new_fac))

        return Polynomial(set(new_terms))

    def __str__(self) -> str:
        return ' + '.join(str(term) for term in self.terms)

    def _validate(self):
        for term in self.terms:
            if isclose(term.factor, 0.0):
                raise ValueError('Polynomial {0} has zero term.'.format(str(self)))


class PolynomialTerm:
    def __init__(self, monomial: Monomial, factor: float):
        self.monomial = monomial
        self.factor = factor

    def __eq__(self, other: 'PolynomialTerm') -> bool:
        return self.monomial == other.monomial and isclose(self.factor, other.factor)

    def __hash__(self) -> int:
        return hash('*puly-term*{0}{1}'.format(str(self.monomial), str(self.factor)))

    def __str__(self) -> str:
        return '{0} * {1}'.format(str(self.factor), str(self.monomial))





