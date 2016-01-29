import pytest

import rxncon.simulation.ode.polynomials as pol


### MONOMIAL ALGEBRA ###
def test_monomial_multiplication_positive_powers():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    z = pol.Symbol('z')

    first_monomial = pol.Monomial({pol.MonomialFactor(x, 2), pol.MonomialFactor(y, 3), pol.MonomialFactor(z, 4)})
    second_monomial = pol.Monomial({pol.MonomialFactor(x, 1), pol.MonomialFactor(y, 1), pol.MonomialFactor(z, 1)})

    assert first_monomial * second_monomial == \
           pol.Monomial({pol.MonomialFactor(x, 3), pol.MonomialFactor(y, 4), pol.MonomialFactor(z, 5)})


def test_monomial_multiplication_mixed_sign_powers():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    z = pol.Symbol('z')

    first_monomial = pol.Monomial({pol.MonomialFactor(x, 2), pol.MonomialFactor(y, 3), pol.MonomialFactor(z, 4)})
    second_monomial = pol.Monomial({pol.MonomialFactor(x, 1), pol.MonomialFactor(y, -3), pol.MonomialFactor(z, 1)})

    assert first_monomial * second_monomial == \
           pol.Monomial({pol.MonomialFactor(x, 3), pol.MonomialFactor(z, 5)})


def test_monomial_multiplication_to_trivial():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    z = pol.Symbol('z')

    first_monomial = pol.Monomial({pol.MonomialFactor(x, 2), pol.MonomialFactor(y, -3), pol.MonomialFactor(z, 4)})
    second_monomial = pol.Monomial({pol.MonomialFactor(x, -2), pol.MonomialFactor(z, -4), pol.MonomialFactor(y, 3)})

    assert first_monomial * second_monomial == pol.Monomial({pol.MonomialFactor(pol.TrivialMulSymbol(), 0)})
    assert second_monomial * first_monomial == pol.TrivialMonomial()


### POLYNOMIAL ALGEBRA ###
def test_polynomial_multiplication(x_plus_y, x_minus_y, x_sq_minus_y_sq, x_sq_plus_two_xy_plus_y_sq):
    assert x_plus_y * x_minus_y == x_sq_minus_y_sq
    assert x_plus_y * x_plus_y == x_sq_plus_two_xy_plus_y_sq


def test_polynomial_addition(x_plus_y, x_minus_y, two_x, two_x_plus_two_y, two_x_minus_two_y):
    assert x_plus_y + x_plus_y == two_x_plus_two_y
    assert x_plus_y + x_minus_y == two_x
    assert x_minus_y + x_minus_y == two_x_minus_two_y


@pytest.fixture
def two_x():
    x = pol.Symbol('x')
    term_two_x = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1)}), 2.0)

    return pol.Polynomial({term_two_x})

@pytest.fixture
def x_plus_y():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    term_x = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1)}), 1.0)
    term_y = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(y, 1)}), 1.0)

    return pol.Polynomial({term_x, term_y})

@pytest.fixture
def two_x_plus_two_y():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    term_two_x = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1)}), 2.0)
    term_two_y = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(y, 1)}), 2.0)

    return pol.Polynomial({term_two_x, term_two_y})

@pytest.fixture
def two_x_minus_two_y():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    term_two_x = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1)}), 2.0)
    term_minus_two_y = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(y, 1)}), -2.0)

    return pol.Polynomial({term_two_x, term_minus_two_y})

@pytest.fixture
def x_minus_y():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    term_x = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1)}), 1.0)
    term_minus_y = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(y, 1)}), -1.0)

    return pol.Polynomial({term_x, term_minus_y})

@pytest.fixture
def x_sq_minus_y_sq():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    term_x_sq = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 2)}), 1.0)
    term_minus_y_sq = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(y, 2)}), -1.0)

    return pol.Polynomial({term_x_sq, term_minus_y_sq})

@pytest.fixture
def x_sq_plus_two_xy_plus_y_sq():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    term_2xy = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1), pol.MonomialFactor(y, 1)}), 2.0)
    term_x_sq = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 2)}), 1.0)
    term_y_sq = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(y, 2)}), 1.0)

    return pol.Polynomial({term_x_sq, term_y_sq, term_2xy})