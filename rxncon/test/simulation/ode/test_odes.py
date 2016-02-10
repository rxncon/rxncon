import pytest
import math

import rxncon.simulation.ode.polynomials as pol
import rxncon.simulation.ode.ode as ode


def test_simple_ode_system(first_lotke_volterra_ode, second_lotke_volterra_ode):
    ode_sys = ode.ODESystem([first_lotke_volterra_ode, second_lotke_volterra_ode])
    assert isinstance(ode_sys, ode.ODESystem)


def test_doubly_defined_time_derivatives(first_lotke_volterra_ode, second_lotke_volterra_ode):
    second_lotke_volterra_ode.time_derivative = pol.Symbol('x')

    with pytest.raises(AssertionError):
        ode_sys = ode.ODESystem([first_lotke_volterra_ode, second_lotke_volterra_ode])


def test_single_ode_to_py_code(first_lotke_volterra_ode, second_lotke_volterra_ode, first_ode_codes, second_ode_codes):
    # Due to commutativity there are four possible code realizations for each ODE.
    assert first_lotke_volterra_ode.time_derivative.name == 'x'
    assert first_lotke_volterra_ode.to_py_code() in first_ode_codes

    assert second_lotke_volterra_ode.time_derivative.name == 'y'
    assert second_lotke_volterra_ode.to_py_code() in second_ode_codes


def test_ode_system_symbol_definition(lotke_volterra_system):
    assert lotke_volterra_system.to_py_code_symbol_defs() == ['x = ys[0]\n', 'y = ys[1]\n']


def test_ode_system_function_definition(lotke_volterra_system):
    assert lotke_volterra_system.to_py_code_return_statement() == 'return [f0, f1]'


def test_ode_odeint_function(lotke_volterra_system):
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


@pytest.fixture
def lotke_volterra_system(first_lotke_volterra_ode, second_lotke_volterra_ode):
    return ode.ODESystem([first_lotke_volterra_ode, second_lotke_volterra_ode])

@pytest.fixture
def first_lotke_volterra_ode():
    # return dx/dt = x - x*y
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    rhs = pol.Polynomial({pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1)}), 1),
                          pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1), pol.MonomialFactor(y, 1)}), -1)})

    return ode.ODE(x, rhs)


@pytest.fixture
def second_lotke_volterra_ode():
    # return dy/dt = x*y - y
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    rhs = pol.Polynomial({pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(y, 1)}), -1),
                          pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1), pol.MonomialFactor(y, 1)}), 1)})

    return ode.ODE(y, rhs)


@pytest.fixture
def first_ode_codes():
    return ['(-1) * (y ** 1 * x ** 1) + (1) * (x ** 1)',  '(1) * (x ** 1) + (-1) * (y ** 1 * x ** 1)',
        '(-1) * (x ** 1 * y ** 1) + (1) * (x ** 1)', '(1) * (x ** 1) + (-1) * (x ** 1 * y ** 1)']


@pytest.fixture
def second_ode_codes():
    return ['(1) * (y ** 1 * x ** 1) + (-1) * (y ** 1)', '(1) * (x ** 1 * y ** 1) + (-1) * (y ** 1)',
        '(-1) * (y ** 1) + (1) * (y ** 1 * x ** 1)', '(-1) * (y ** 1) + (1) * (x ** 1 * y ** 1)']