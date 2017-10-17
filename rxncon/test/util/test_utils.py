from rxncon.util.utils import elems_eq


def test_elems_eq_not_nested() -> None:
    assert elems_eq([1, 2, 3], [3, 1, 2])
    assert not elems_eq([1, 2], [1, 2, 3])

    assert elems_eq([], [])       # type: ignore
    assert not elems_eq([1], [])  # type: ignore


def test_elems_eq_nested() -> None:
    assert elems_eq([[1, 2, 3], [4, 5]], [[2, 1, 3], [5, 4]])
    assert elems_eq([[4, 5], [3, 2, 1]], [[2, 1, 3], [5, 4]])

    assert elems_eq([[], [1, 2, 3]], [[3, 2, 1], []])  # type: ignore
    assert elems_eq([[], []], [[], []])                # type: ignore
