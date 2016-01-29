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
