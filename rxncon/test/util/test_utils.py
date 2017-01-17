from rxncon.util.utils import elems_eq, current_function_name


def test_elems_eq_not_nested():
    assert elems_eq([1, 2, 3], [3, 1, 2])
    assert not elems_eq([1, 2], [1, 2, 3])

    assert elems_eq([], [])
    assert not elems_eq([1], [])


def test_elems_eq_nested():
    assert elems_eq([[1, 2, 3], [4, 5]], [[2, 1, 3], [5, 4]])
    assert elems_eq([[4, 5], [3, 2, 1]], [[2, 1, 3], [5, 4]])

    assert elems_eq([[], [1, 2, 3]], [[3, 2, 1], []])
    assert elems_eq([[], []], [[], []])


def test_current_function_name():
    assert current_function_name(colored=False) == 'test_current_function_name'
