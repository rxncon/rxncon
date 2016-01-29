from typing import List, Callable, Set
from collections import defaultdict

import typecheck as tc


TRIVIAL_SYMBOL = '1'


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
        return '*symb*{0}'.format(self.name)

    @property
    def is_trivial(self):
        return self.name == TRIVIAL_SYMBOL


def TrivialSymbol():
    return Symbol(TRIVIAL_SYMBOL)


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
            return Monomial({MonomialFactor(TrivialSymbol(), 0)})

        return Monomial(set(new_factors))


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
        if self.power == 0 and not self.symbol.is_trivial:
            raise ValueError('Only TrivialSymbol may be raised to the zero-th power.')

        elif self.symbol.is_trivial and not self.power == 0:
            raise ValueError('TrivialSymbol may only be raised to the zero-th power.')







