import pytest
import rxncon.core.specification as spec


@pytest.fixture
def specifications_input():
    return [[spec.ProteinSpecification("C", None, None, None), spec.ProteinSpecification("B", None, None, None), spec.ProteinSpecification("A", None, None, None)],
            [spec.ProteinSpecification("A", "dC", None, None), spec.ProteinSpecification("A", "dB", None, None), spec.ProteinSpecification("A", "dA", None, None)],
            [spec.ProteinSpecification("A", "d", "sC", None), spec.ProteinSpecification("A", "d", "sB", None), spec.ProteinSpecification("A", "d", "sA", None)],
            [spec.ProteinSpecification("A", None, None, "rC"), spec.ProteinSpecification("A", None, None, "rB"), spec.ProteinSpecification("A", None, None, "rA")],
            [spec.ProteinSpecification("A", "dC", None, None), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", "dA", None, None)],
            [spec.ProteinSpecification("A", "d", "sC", None), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", "d", "sA", None)],
            [spec.ProteinSpecification("A", None, None, "rC"), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, "rA")],
            [spec.ProteinSpecification("C", "dC", "sC", "rC"), spec.ProteinSpecification("B", "dB", "sB", "rB"), spec.ProteinSpecification("A", "dA", "sA", "rA")],
            [spec.ProteinSpecification("A", None, None, "rC"), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, None)],

            [spec.RnaSpecification("C", None, None, None), spec.RnaSpecification("B", None, None, None), spec.RnaSpecification("A", None, None, None)],
            [spec.RnaSpecification("A", "dC", None, None), spec.RnaSpecification("A", "dB", None, None), spec.RnaSpecification("A", "dA", None, None)],
            [spec.RnaSpecification("A", "d", "sC", None), spec.RnaSpecification("A", "d", "sB", None), spec.RnaSpecification("A", "d", "sA", None)],
            [spec.RnaSpecification("A", None, None, "rC"), spec.RnaSpecification("A", None, None, "rB"), spec.RnaSpecification("A", None, None, "rA")],
            [spec.RnaSpecification("A", "dC", None, None), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", "dA", None, None)],
            [spec.RnaSpecification("A", "d", "sC", None), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", "d", "sA", None)],
            [spec.RnaSpecification("A", None, None, "rC"), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, "rA")],
            [spec.RnaSpecification("C", "dC", "sC", "rC"), spec.RnaSpecification("B", "dB", "sB", "rB"), spec.RnaSpecification("A", "dA", "sA", "rA")],
            [spec.RnaSpecification("A", None, None, "rC"), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None)]
            ]


@pytest.fixture
def expected_specification_ordering():
    return [[spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("B", None, None, None), spec.ProteinSpecification("C", None, None, None)],
            [spec.ProteinSpecification("A", "dA", None, None), spec.ProteinSpecification("A", "dB", None, None), spec.ProteinSpecification("A", "dC", None, None)],
            [spec.ProteinSpecification("A", "d", "sA", None), spec.ProteinSpecification("A", "d", "sB", None), spec.ProteinSpecification("A", "d", "sC", None)],
            [spec.ProteinSpecification("A", None, None, "rA"), spec.ProteinSpecification("A", None, None, "rB"), spec.ProteinSpecification("A", None, None, "rC")],
            [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", "dA", None, None), spec.ProteinSpecification("A", "dC", None, None)],
            [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", "d", "sA", None), spec.ProteinSpecification("A", "d", "sC", None)],
            [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, "rA"), spec.ProteinSpecification("A", None, None, "rC")],
            [spec.ProteinSpecification("A", "dA", "sA", "rA"), spec.ProteinSpecification("B", "dB", "sB", "rB"), spec.ProteinSpecification("C", "dC", "sC", "rC")],
            [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, "rC")],

            [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("B", None, None, None), spec.RnaSpecification("C", None, None, None)],
            [spec.RnaSpecification("A", "dA", None, None), spec.RnaSpecification("A", "dB", None, None), spec.RnaSpecification("A", "dC", None, None)],
            [spec.RnaSpecification("A", "d", "sA", None), spec.RnaSpecification("A", "d", "sB", None), spec.RnaSpecification("A", "d", "sC", None)],
            [spec.RnaSpecification("A", None, None, "rA"), spec.RnaSpecification("A", None, None, "rB"), spec.RnaSpecification("A", None, None, "rC")],
            [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", "dA", None, None), spec.RnaSpecification("A", "dC", None, None)],
            [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", "d", "sA", None), spec.RnaSpecification("A", "d", "sC", None)],
            [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, "rA"), spec.RnaSpecification("A", None, None, "rC")],
            [spec.RnaSpecification("A", "dA", "sA", "rA"), spec.RnaSpecification("B", "dB", "sB", "rB"), spec.RnaSpecification("C", "dC", "sC", "rC")],
            [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, "rC")]
           ]

@pytest.fixture
def specification_classes():
    return [[spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None)],
            [spec.ProteinSpecification("A", "dA", None, None), spec.RnaSpecification("A", "dB", None, None)],
            [spec.ProteinSpecification("A", "d", "sA", None), spec.RnaSpecification("A", "d", "sB", None)],
            [spec.ProteinSpecification("A", None, None, "rA"), spec.RnaSpecification("A", None, None, "rB")],
            [spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", "dA", None, None)],
            [spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", "d", "sA", None)],
            [spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, "rA")],
            [spec.ProteinSpecification("A", "dA", "sA", "rA"), spec.RnaSpecification("B", "dB", "sB", "rB")],
            [spec.RnaSpecification("B", None, None, None), spec.ProteinSpecification("B", None, None, None)],
            [spec.ProteinSpecification("A", "dA", "sA", "rA"), spec.RnaSpecification("B", "dB", "sB", "rB"), spec.RnaSpecification("A", "dB", "sB", None)]
            ]


@pytest.fixture
def expected_specification_class_sorting():
    return [[spec.RnaSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, None)],
            [spec.RnaSpecification("A", "dB", None, None), spec.ProteinSpecification("A", "dA", None, None)],
            [spec.RnaSpecification("A", "d", "sB", None), spec.ProteinSpecification("A", "d", "sA", None)],
            [spec.RnaSpecification("A", None, None, "rB"), spec.ProteinSpecification("A", None, None, "rA")],
            [spec.RnaSpecification("A", "dA", None, None), spec.ProteinSpecification("A", None, None, None)],
            [spec.RnaSpecification("A", "d", "sA", None), spec.ProteinSpecification("A", None, None, None)],
            [spec.RnaSpecification("A", None, None, "rA"), spec.ProteinSpecification("A", None, None, None)],
            [spec.RnaSpecification("B", "dB", "sB", "rB"), spec.ProteinSpecification("A", "dA", "sA", "rA")],
            [spec.RnaSpecification("B", None, None, None), spec.ProteinSpecification("B", None, None, None)],
            [spec.RnaSpecification("A", "dB", "sB", None), spec.RnaSpecification("B", "dB", "sB", "rB"), spec.ProteinSpecification("A", "dA", "sA", "rA")]
            ]


def test_specification_internal_sorting(specifications_input, expected_specification_ordering):
    for i, specifications in enumerate(specifications_input):
        assert sorted(specifications) == expected_specification_ordering[i]


def test_specification_sorting(specification_classes, expected_specification_class_sorting):
    for i, specifications in enumerate(specification_classes):
        assert sorted(specifications) == expected_specification_class_sorting[i]

