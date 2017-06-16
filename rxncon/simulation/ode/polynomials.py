from typing import List, Set, Union, Dict
from collections import defaultdict
from math import isclose


TRIVIAL_MUL_SYMBOL = '1'
TRIVIAL_ADD_SYMBOL = '0'


### SYMBOLS, i.e. x, y, z ###
class Symbol:
    def __init__(self, name: str) -> None:
        self.name = name

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Symbol):
            return NotImplemented
        return self.name == other.name

    def __hash__(self) -> int:
        return hash(self.name)

    def __str__(self) -> str:
        return '{0}'.format(self.name)

    @property
    def is_trivial_mul(self) -> bool:
        return self.name == TRIVIAL_MUL_SYMBOL

    @property
    def is_trivial_add(self) -> bool:
        return self.name == TRIVIAL_ADD_SYMBOL


def TrivialMulSymbol() -> Symbol:
    return Symbol(TRIVIAL_MUL_SYMBOL)


def TrivialAddSymbol() -> Symbol:
    return Symbol(TRIVIAL_ADD_SYMBOL)


### MONOMIALS, i.e. (x^3)*(y^4) ###
class Monomial:
    def __init__(self, factors: Set['MonomialFactor']) -> None:
        self.factors = factors

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Monomial):
            return NotImplemented
        return self.factors == other.factors

    def __hash__(self) -> int:
        return hash('*monomial*{0}'.join(str(x) for x in self.factors))

    def __mul__(self, other: 'Monomial') -> 'Monomial':
        sym_to_power = defaultdict(int)  # type: Dict[Symbol, int]

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

    @property
    def is_constant(self) -> bool:
        return self == TrivialMonomial()

    @property
    def symbols(self) -> Set[Symbol]:
        return {x.symbol for x in self.factors}


class MonomialFactor:
    def __init__(self, symbol: Symbol, power: int) -> None:
        self.symbol = symbol
        self.power = power
        self._validate()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, MonomialFactor):
            return NotImplemented
        return self.symbol == other.symbol and self.power == other.power

    def __hash__(self) -> int:
        return hash('*mono-fac*{0}^{1}'.format(str(self.symbol), str(self.power)))

    def __str__(self) -> str:
        if self.power < 0:
            return '{0}^({1})'.format(self.symbol, self.power)
        elif self.power == 0:
            return '1'
        else:
            return '{0}^{1}'.format(self.symbol, self.power)

    def _validate(self) -> None:
        if self.power == 0 and not self.symbol.is_trivial_mul:
            raise ValueError('Only TrivialSymbol may be raised to the zero-th power.')
        elif self.symbol.is_trivial_mul and not self.power == 0:
            raise ValueError('TrivialSymbol may only be raised to the zero-th power.')


def TrivialMonomial() -> Monomial:
    return Monomial({MonomialFactor(TrivialMulSymbol(), 0)})


### POLYNOMIALS, i.e. 3*(x^2)*(y^4) - 2*(z^8) ###
class Polynomial:
    def __init__(self, terms: Set['PolynomialTerm']) -> None:
        self.terms = terms
        self._validate()

    def __add__(self, other: Union['Polynomial', float]) -> 'Polynomial':
        new_monomial_to_factor = defaultdict(lambda: 0.0)  # type: Dict[Monomial, float]

        if isinstance(other, Polynomial):
            for term in self.terms:
                new_monomial_to_factor[term.monomial] = term.factor
            for term in other.terms:
                new_monomial_to_factor[term.monomial] += term.factor

        elif isinstance(other, int) or isinstance(other, float):
            for term in self.terms:
                new_monomial_to_factor[term.monomial] = term.factor
            new_monomial_to_factor[TrivialMonomial()] += other

        new_terms = []

        for new_mon, new_fac in new_monomial_to_factor.items():
            if not isclose(new_fac, 0.0):
                new_terms.append(PolynomialTerm(new_mon, new_fac))

        return Polynomial(set(new_terms))

    def __sub__(self, other: Union['Polynomial', float]) -> 'Polynomial':
        if isinstance(other, Polynomial):
            return self + (-1 * other)
        elif isinstance(other, float) or isinstance(other, int):
            return self + (-1 * other)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Polynomial):
            return NotImplemented
        return self.terms == other.terms

    def __mul__(self, mult_by: Union['Polynomial', float]) -> 'Polynomial':
        if isinstance(mult_by, Polynomial):
            this_monomial_to_factor = {term.monomial: term.factor for term in self.terms}
            that_monomial_to_factor = {term.monomial: term.factor for term in mult_by.terms}

            new_monomial_to_factor = defaultdict(lambda: 0.0)  # type: Dict[Monomial, float]

            for this_mon, this_fac in this_monomial_to_factor.items():
                for that_mon, that_fac in that_monomial_to_factor.items():
                    new_monomial_to_factor[this_mon * that_mon] += this_fac * that_fac

            new_terms = []

            for new_mon, new_fac in new_monomial_to_factor.items():
                if not isclose(new_fac, 0.0):
                    new_terms.append(PolynomialTerm(new_mon, new_fac))

            return Polynomial(set(new_terms))
        elif isinstance(mult_by, float) or isinstance(mult_by, int):
            new_terms = []

            for term in self.terms:
                new_terms.append(PolynomialTerm(term.monomial, term.factor * mult_by))

            return Polynomial(set(new_terms))

    def __rmul__(self, other: Union['Polynomial', float]) -> 'Polynomial':  # type: ignore
        return self * other

    def __str__(self) -> str:
        return ' + '.join(str(term) for term in self.terms)

    @property
    def symbols(self) -> Set[Symbol]:
        syms = []  # type: List[Symbol]
        for term in self.terms:
            syms += list(term.symbols)

        return set(syms)

    def _validate(self) -> None:
        for term in self.terms:
            if isclose(term.factor, 0.0):
                raise ValueError('Polynomial {0} has zero term.'.format(str(self)))


class PolynomialTerm:
    def __init__(self, monomial: Monomial, factor: Union[float]) -> None:
        self.monomial = monomial
        self.factor = factor

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PolynomialTerm):
            return NotImplemented
        return self.monomial == other.monomial and isclose(self.factor, other.factor)

    def __hash__(self) -> int:
        return hash('*puly-term*{0}{1}'.format(str(self.monomial), str(self.factor)))

    def __str__(self) -> str:
        return '{0} * {1}'.format(str(self.factor), str(self.monomial))

    @property
    def is_constant(self) -> bool:
        return self.monomial.is_constant

    @property
    def symbols(self) -> Set[Symbol]:
        return self.monomial.symbols
