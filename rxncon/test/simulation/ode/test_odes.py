import pytest

import rxncon.simulation.ode.polynomials as pol
import rxncon.simulation.ode.ode as ode


def test_simple_ode_system(first_lotke_volterra_ode, second_lotke_volterra_ode):
    ode_sys = ode.ODESystem([first_lotke_volterra_ode, second_lotke_volterra_ode])
    assert isinstance(ode_sys, ode.ODESystem)


def test_doubly_defined_time_derivatives(first_lotke_volterra_ode, second_lotke_volterra_ode):
    second_lotke_volterra_ode.time_derivative = pol.Symbol('x')

    with pytest.raises(AssertionError):
        ode_sys = ode.ODESystem([first_lotke_volterra_ode, second_lotke_volterra_ode])


def test_single_ode_to_py_code(first_lotke_volterra_ode, second_lotke_volterra_ode):
    # The order of the terms in each ODE is immaterial, so there are two possible code realizations for the ODEs.
    assert first_lotke_volterra_ode.to_py_code() == 'd_t_x = (-1) * (y ** 1 * x ** 1) + (1) * (x ** 1)' or \
        first_lotke_volterra_ode.to_py_code() == 'd_t_x = (1) * (x ** 1) + (-1) * (y ** 1 * x ** 1)'
    assert second_lotke_volterra_ode.to_py_code() == 'd_t_y = (1) * (y ** 1 * x ** 1) + (-1) * (y ** 1)' or \
        second_lotke_volterra_ode.to_py_code() == 'd_t_y = (-1) * (y ** 1) + (1) * (y ** 1 * x ** 1)'

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
