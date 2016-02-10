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


def test_monomial_symbols():
    x = pol.Symbol('x')
    y = pol.Symbol('y')
    z = pol.Symbol('z')

    pol.Monomial({pol.MonomialFactor(x, 2), pol.MonomialFactor(y, 3), pol.MonomialFactor(z, 4)})

    assert pol.Monomial({pol.MonomialFactor(x, 2), pol.MonomialFactor(y, 3), pol.MonomialFactor(z, 4)}).symbols == \
        {x, y, z}

    assert pol.Monomial({pol.MonomialFactor(x, 2), pol.MonomialFactor(y, 3)}).symbols == {x, y}


### POLYNOMIAL ALGEBRA ###
def test_polynomial_multiplication_by_polynomial(x_plus_y, x_minus_y, x_sq_minus_y_sq, x_sq_plus_two_xy_plus_y_sq):
    assert x_plus_y * x_minus_y == x_sq_minus_y_sq
    assert x_plus_y * x_plus_y == x_sq_plus_two_xy_plus_y_sq


def test_polynomial_multiplication_by_int(x_plus_y, two_x_plus_two_y):
    assert x_plus_y * 2 == two_x_plus_two_y
    assert 2 * x_plus_y == two_x_plus_two_y


def test_polynomial_addition(x_plus_y, x_minus_y, two_x, two_x_plus_two_y, two_x_minus_two_y):
    assert x_plus_y + x_plus_y == two_x_plus_two_y
    assert x_plus_y + x_minus_y == two_x
    assert x_minus_y + x_minus_y == two_x_minus_two_y


def test_polynomial_addition_by_scalar(one_x, two_x, const_one):
    assert 2 * (one_x + 1) == 2 * (one_x + const_one) == two_x + 2 == two_x + 2 * const_one == two_x + const_one * 2


def test_polynomial_subtraction(one_x, one_y, two_x_minus_two_y, const_one):
    assert 2 * one_x - 2 * one_y == two_x_minus_two_y


def test_polynomial_symbols(two_x, x_sq_minus_y_sq):
    x = pol.Symbol('x')
    y = pol.Symbol('y')

    assert two_x.symbols == {x}
    assert x_sq_minus_y_sq.symbols == {x, y}


### TEST FIXTURES ###
@pytest.fixture
def const_one():
    x = pol.TrivialMulSymbol()
    one = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 0)}), 1.0)

    return pol.Polynomial({one})

@pytest.fixture
def one_x():
    x = pol.Symbol('x')
    one_x = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(x, 1)}), 1.0)

    return pol.Polynomial({one_x})


@pytest.fixture
def one_y():
    y = pol.Symbol('y')
    one_y = pol.PolynomialTerm(pol.Monomial({pol.MonomialFactor(y, 1)}), 1.0)

    return pol.Polynomial({one_y})


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