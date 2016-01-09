from abc import ABCMeta

import rxncon.core.state as sta


class Effector:
    __metaclass__ = ABCMeta


class StateEffector(Effector):
    def __init__(self, expr: sta.State):
        assert isinstance(expr, sta.State)
        self.expr = expr

    def __str__(self) -> str:
        return 'StateEffector({})'.format(self.expr)

    def __eq__(self, other: Effector) -> bool:
        assert isinstance(other, Effector)

        if isinstance(other, StateEffector):
            return self.expr == other.expr

        else:
            return False


class NotEffector(Effector):
    def __init__(self, expr: Effector):
        assert isinstance(expr, Effector)
        self.expr = expr

    def __str__(self) -> str:
        return 'NotEffector({})'.format(self.expr)

    def __eq__(self, other: Effector) -> bool:
        assert isinstance(other, Effector)

        if isinstance(other, NotEffector):
            return self.expr == other.expr

        else:
            return False


class BinaryEffector(Effector):
    __metaclass__ = ABCMeta

    def __init__(self, left_expr: Effector, right_expr: Effector):
        self.left_expr = left_expr
        self.right_expr = right_expr

    def __eq__(self, other: Effector) -> bool:
        assert isinstance(other, Effector)

        if isinstance(other, BinaryEffector):
            return self.left_expr.__class__ == other.left_expr.__class__ and self.right_expr.__class__ == other.right_expr.__class__ and \
                self.left_expr == other.left_expr and self.right_expr == other.right_expr

        else:
            return False


class AndEffector(BinaryEffector):
    def __str__(self) -> str:
        return 'AndEffector({0}, {1})'.format(self.left_expr, self.right_expr)


class OrEffector(BinaryEffector):
    def __str__(self) -> str:
        return 'OrEffector({0}, {1})'.format(self.left_expr, self.right_expr)