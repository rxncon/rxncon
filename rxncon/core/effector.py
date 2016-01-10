from abc import ABCMeta

import rxncon.core.state as sta


class Effector:
    """Effector is a tree data structure with as its leafs StateEffector objects, which hold State objects. The non-leaf
    node types NotEffector, AndEffector and OrEffector represent boolean NOT, AND and OR."""
    __metaclass__ = ABCMeta

    @property
    def name(self):
        """The optional name of the Effector. When an Effector tree is constructed through the read-in and parsing of
        ContingencyListEntry objects, the name of the boolean Effector objects is set to the value with which they appear
        in the contingency list."""
        try:
            return self._name

        except AttributeError:
            return None

    @name.setter
    def name(self, value: str):
        assert isinstance(value, str)
        self._name = value


class StateEffector(Effector):
    def __init__(self, expr: sta.State):
        assert isinstance(expr, sta.State)
        self.expr = expr

    def __str__(self) -> str:
        return 'StateEffector({})'.format(self.expr)

    def __eq__(self, other: Effector) -> bool:
        assert isinstance(other, Effector)

        if isinstance(other, StateEffector):
            return self.expr == other.expr and self.name == other.name

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
            return self.expr == other.expr and self.name == other.name

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
                self.left_expr == other.left_expr and self.right_expr == other.right_expr and self.left_expr.name == other.left_expr.name and \
                self.right_expr.name == other.right_expr.name

        else:
            return False


class AndEffector(BinaryEffector):
    def __str__(self) -> str:
        return 'AndEffector({0}, {1})'.format(self.left_expr, self.right_expr)


class OrEffector(BinaryEffector):
    def __str__(self) -> str:
        return 'OrEffector({0}, {1})'.format(self.left_expr, self.right_expr)