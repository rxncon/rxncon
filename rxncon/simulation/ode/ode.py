import typecheck as tc

import rxncon.simulation.ode.polynomials as pol


class ODE:
    @tc.typecheck
    def __init__(self, time_derivative: pol.Symbol, polynomial: pol.Polynomial):
        self.time_derivative = time_derivative
        self.polynomial = polynomial

    @property
    def symbols(self):
        return self.polynomial.symbols.union({self.time_derivative})

