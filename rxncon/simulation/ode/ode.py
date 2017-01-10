from typing import List, Set as Set, Callable

from rxncon.simulation.ode.polynomials import Polynomial, PolynomialTerm, Symbol, Monomial, MonomialFactor


PythonStatement = str


class ODE:
    def __init__(self, time_derivative: Symbol, polynomial: Polynomial) -> None:
        self.time_derivative = time_derivative
        self.polynomial = polynomial

    @property
    def symbols(self) -> Set[Symbol]:
        return self.polynomial.symbols.union({self.time_derivative})

    def to_py_code(self) -> PythonStatement:
        rhs = ' + '.join(_polynomial_term_to_py_code(x) for x in self.polynomial.terms)

        return '{0}'.format(rhs)


class ODESystem:
    name_func = 'f'
    name_y = 'ys'
    name_t = 't'

    def __init__(self, odes: List[ODE]) -> None:
        self.odes = odes
        self._validate()

    @property
    def odeint_function(self) -> Callable:
        indent = '    '
        body_statements = self.to_py_code_symbol_defs() + self.to_py_code_function_defs() + [self.to_py_code_return_statement()]
        body_string = indent + indent.join(body_statements)

        exec('def {0}({1}, {2}):\n{3}'.format(self.name_func, self.name_y, self.name_t, body_string))

        return eval(self.name_func)

    def to_py_code_symbol_defs(self) -> List[PythonStatement]:
        symbol_defs = []

        for i, ode in enumerate(self.odes):
            symbol_defs.append('{0} = {1}[{2}]\n'.format(ode.time_derivative.name, self.name_y, i))

        return symbol_defs

    def to_py_code_function_defs(self) -> List[PythonStatement]:
        function_defs = []

        for i, ode in enumerate(self.odes):
            function_defs.append('{0}{1} = {2}\n'.format(self.name_func, i, ode.to_py_code()))

        return function_defs

    def to_py_code_return_statement(self) -> PythonStatement:
        func_names = ['{0}{1}'.format(self.name_func, i) for i, ode in enumerate(self.odes)]
        func_string = ', '.join(func_names)

        return 'return [{0}]'.format(func_string)

    def _validate(self) -> None:
        if not len(set(x.time_derivative for x in self.odes)) == len(self.odes):
            raise AssertionError('Multiply defined time_derivative terms in ODESystem {}'.format(self.odes))


def _polynomial_term_to_py_code(polterm: PolynomialTerm) -> PythonStatement:
    if polterm.monomial.is_constant:
        return '{0}'.format(polterm.factor)
    else:
        return '({0}) * ({1})'.format(polterm.factor, _monomial_to_py_code(polterm.monomial))


def _monomial_to_py_code(monomial: Monomial) -> PythonStatement:
    return ' * '.join(_monomial_factor_to_py_code(x) for x in monomial.factors)


def _monomial_factor_to_py_code(monomial_factor: MonomialFactor) -> PythonStatement:
    return '{0} ** {1}'.format(monomial_factor.symbol.name, monomial_factor.power)