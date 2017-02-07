import pytest
import math
from typing import List

from rxncon.simulation.ode.polynomials import Polynomial, PolynomialTerm, Symbol, Monomial, MonomialFactor
from rxncon.simulation.ode.ode import ODE, ODESystem


def test_simple_ode_system(first_lotke_volterra_ode: ODE, second_lotke_volterra_ode: ODE) -> None:
    ode_sys = ODESystem([first_lotke_volterra_ode, second_lotke_volterra_ode])
    assert isinstance(ode_sys, ODESystem)


def test_doubly_defined_time_derivatives(first_lotke_volterra_ode: ODE, second_lotke_volterra_ode: ODE) -> None:
    second_lotke_volterra_ode.time_derivative = Symbol('x')

    with pytest.raises(AssertionError):
        ODESystem([first_lotke_volterra_ode, second_lotke_volterra_ode])


def test_single_ode_to_py_code(first_lotke_volterra_ode: ODE, second_lotke_volterra_ode: ODE,
                               first_ode_codes: List[str], second_ode_codes: List[str]) -> None:
    # Due to commutativity there are four possible code realizations for each ODE.
    assert first_lotke_volterra_ode.time_derivative.name == 'x'
    assert first_lotke_volterra_ode.to_py_code() in first_ode_codes

    assert second_lotke_volterra_ode.time_derivative.name == 'y'
    assert second_lotke_volterra_ode.to_py_code() in second_ode_codes


def test_ode_system_symbol_definition(lotke_volterra_system: ODESystem) -> None:
    assert lotke_volterra_system.to_py_code_symbol_defs() == ['x = ys[0]\n', 'y = ys[1]\n']


def test_ode_system_function_definition(lotke_volterra_system: ODESystem) -> None:
    assert lotke_volterra_system.to_py_code_return_statement() == 'return [f0, f1]'


def test_ode_odeint_function(lotke_volterra_system: ODESystem) -> None:
    f = lotke_volterra_system.odeint_function

    x_y_to_expected_f_x_f_y = [
        ([0.2, 0.3], [0.14, -0.24]),
        ([0.0, 0.0], [0.0, 0.0]),
        ([0.1, 0.1], [0.09, -0.09]),
        ([0.5, 1.0], [0.0, -0.5]),
        ([1.0, 0.5], [0.5, 0.0]),
        ([2.0, 0.1], [1.8, 0.1]),
        ([-3.0, 2.0], [3.0, -8.0]),
        ([0.3, -0.8], [0.54, 0.56])
    ]

    # Time independent, so t == t0
    t0 = 0.0
    for xs, fs in x_y_to_expected_f_x_f_y:
        assert math.isclose(f(xs, t0)[0], fs[0]) and math.isclose(f(xs, t0)[1], fs[1])


def test_ode_odeint(lotke_volterra_system: ODESystem) -> None:
    import numpy as np
    import scipy.integrate as integ

    f = lotke_volterra_system.odeint_function
    y0 = [4.0, 5.0]
    t = np.linspace(0.0, 15.0, 100)

    solution = integ.odeint(f, y0, t)

    # Compare this solution with one generated in Mathematica:
    expected_index_to_x_to_y = {
        0:  [4.0, 5.0],
        20: [0.00796014, 0.532128],
        40: [0.0986716, 0.0284035],
        60: [1.97129, 0.00907175],
        80: [0.017887, 3.09127],
        99: [0.016797, 0.178667]
    }

    # The tolerance is quite high, since the numbers coming from Mathematica don't have such a high precision.
    for index, x_y in expected_index_to_x_to_y.items():
        assert math.isclose(solution[index, 0], x_y[0], abs_tol=1E-5) and \
               math.isclose(solution[index, 1], x_y[1], abs_tol=1E-5)


@pytest.fixture
def lotke_volterra_system(first_lotke_volterra_ode: ODE, second_lotke_volterra_ode: ODE) -> ODESystem:
    return ODESystem([first_lotke_volterra_ode, second_lotke_volterra_ode])

@pytest.fixture
def first_lotke_volterra_ode() -> ODE:
    # return dx/dt = x - x*y
    x = Symbol('x')
    y = Symbol('y')
    rhs = Polynomial({PolynomialTerm(Monomial({MonomialFactor(x, 1)}), 1),
                      PolynomialTerm(Monomial({MonomialFactor(x, 1), MonomialFactor(y, 1)}), -1)})

    return ODE(x, rhs)


@pytest.fixture
def second_lotke_volterra_ode() -> ODE:
    # return dy/dt = x*y - y
    x = Symbol('x')
    y = Symbol('y')
    rhs = Polynomial({PolynomialTerm(Monomial({MonomialFactor(y, 1)}), -1),
                      PolynomialTerm(Monomial({MonomialFactor(x, 1), MonomialFactor(y, 1)}), 1)})

    return ODE(y, rhs)


@pytest.fixture
def first_ode_codes() -> List[str]:
    return ['(-1) * (y ** 1 * x ** 1) + (1) * (x ** 1)',  '(1) * (x ** 1) + (-1) * (y ** 1 * x ** 1)',
            '(-1) * (x ** 1 * y ** 1) + (1) * (x ** 1)', '(1) * (x ** 1) + (-1) * (x ** 1 * y ** 1)']


@pytest.fixture
def second_ode_codes() -> List[str]:
    return ['(1) * (y ** 1 * x ** 1) + (-1) * (y ** 1)', '(1) * (x ** 1 * y ** 1) + (-1) * (y ** 1)',
            '(-1) * (y ** 1) + (1) * (y ** 1 * x ** 1)', '(-1) * (y ** 1) + (1) * (x ** 1 * y ** 1)']