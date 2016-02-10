import typecheck as tc
from typing import List

import rxncon.simulation.ode.polynomials as pol



class ODE:
    prefix = 'd_t_'

    @tc.typecheck
    def __init__(self, time_derivative: pol.Symbol, polynomial: pol.Polynomial):
        self.time_derivative = time_derivative
        self.polynomial = polynomial

    @property
    def symbols(self):
        return self.polynomial.symbols.union({self.time_derivative})

    def to_py_code(self):
        lhs = self.prefix + self.time_derivative.name
        rhs = ' + '.join(_polynomial_term_to_py_code(x) for x in self.polynomial.terms)

        return '{0} = {1}'.format(lhs, rhs)


class ODESystem:
    def __init__(self, odes: List[ODE]):
        self.odes = odes
        self._validate()

    def _validate(self):
        if not len(set(x.time_derivative for x in self.odes)) == len(self.odes):
            raise AssertionError('Multiply defined time_derivative terms in ODESystem {}'.format(self.odes))


def _polynomial_term_to_py_code(polterm: pol.PolynomialTerm):
    if polterm.monomial.is_constant:
        return '{0}'.format(polterm.factor)
    else:
        return '({0}) * ({1})'.format(polterm.factor, _monomial_to_py_code(polterm.monomial))


def _monomial_to_py_code(monomial: pol.Monomial):
    return ' * '.join(_monomial_factor_to_py_code(x) for x in monomial.factors)


def _monomial_factor_to_py_code(monomial_factor: pol.MonomialFactor):
    return '{0} ** {1}'.format(monomial_factor.symbol.name, monomial_factor.power)