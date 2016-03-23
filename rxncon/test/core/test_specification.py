import pytest
import rxncon.core.specification as spec


@pytest.fixture
def specifications_input():
    return [[spec.Specification("C", None, None, None), spec.Specification("B", None, None, None), spec.Specification("A", None, None, None)],
            [spec.Specification("A", "dC", None, None), spec.Specification("A", "dB", None, None), spec.Specification("A", "dA", None, None)],
            [spec.Specification("A", "d", "sC", None), spec.Specification("A", "d", "sB", None), spec.Specification("A", "d", "sA", None)],
            [spec.Specification("A", None, None, "rC"), spec.Specification("A", None, None, "rB"), spec.Specification("A", None, None, "rA")],
            [spec.Specification("A", "dC", None, None), spec.Specification("A", None, None, None), spec.Specification("A", "dA", None, None)],
            [spec.Specification("A", "d", "sC", None), spec.Specification("A", None, None, None), spec.Specification("A", "d", "sA", None)],
            [spec.Specification("A", None, None, "rC"), spec.Specification("A", None, None, None), spec.Specification("A", None, None, "rA")],
            [spec.Specification("C", "dC", "sC", "rC"), spec.Specification("B", "dB", "sB", "rB"), spec.Specification("A", "dA", "sA", "rA")]
            ]


@pytest.fixture
def expected_specification_ordering():
    return [[spec.Specification("A", None, None, None), spec.Specification("B", None, None, None), spec.Specification("C", None, None, None)],
            [spec.Specification("A", "dA", None, None), spec.Specification("A", "dB", None, None), spec.Specification("A", "dC", None, None)],
            [spec.Specification("A", "d", "sA", None), spec.Specification("A", "d", "sB", None), spec.Specification("A", "d", "sC", None)],
            [spec.Specification("A", None, None, "rA"), spec.Specification("A", None, None, "rB"), spec.Specification("A", None, None, "rC")],
            [spec.Specification("A", None, None, None), spec.Specification("A", "dA", None, None), spec.Specification("A", "dC", None, None)],
            [spec.Specification("A", None, None, None), spec.Specification("A", "d", "sA", None), spec.Specification("A", "d", "sC", None)],
            [spec.Specification("A", None, None, None), spec.Specification("A", None, None, "rA"), spec.Specification("A", None, None, "rC")],
            [spec.Specification("A", "dA", "sA", "rA"), spec.Specification("B", "dB", "sB", "rB"), spec.Specification("C", "dC", "sC", "rC")]
           ]


def test_specification_sorting(specifications_input, expected_specification_ordering):
    for i, specifications in enumerate(specifications_input):
        assert sorted(specifications) == expected_specification_ordering[i]