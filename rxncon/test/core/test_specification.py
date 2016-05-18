import pytest
from collections import namedtuple
import rxncon.core.specification as spec

OrderTestCase = namedtuple('OrderTestCase', ['to_sort', 'expected_order'])
@pytest.fixture
def the_case_specifications_ordering():
    return [
            OrderTestCase([spec.ProteinSpecification("C", None, None, None), spec.ProteinSpecification("B", None, None, None), spec.ProteinSpecification("A", None, None, None)],
                          [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("B", None, None, None), spec.ProteinSpecification("C", None, None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", "dC", None, None), spec.ProteinSpecification("A", "dB", None, None), spec.ProteinSpecification("A", "dA", None, None)],
                          [spec.ProteinSpecification("A", "dA", None, None), spec.ProteinSpecification("A", "dB", None, None), spec.ProteinSpecification("A", "dC", None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", "d", "sC", None), spec.ProteinSpecification("A", "d", "sB", None), spec.ProteinSpecification("A", "d", "sA", None)],
                          [spec.ProteinSpecification("A", "d", "sA", None), spec.ProteinSpecification("A", "d", "sB", None), spec.ProteinSpecification("A", "d", "sC", None)]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, "rC"), spec.ProteinSpecification("A", None, None, "rB"), spec.ProteinSpecification("A", None, None, "rA")],
                          [spec.ProteinSpecification("A", None, None, "rA"), spec.ProteinSpecification("A", None, None, "rB"), spec.ProteinSpecification("A", None, None, "rC")]),


            OrderTestCase([spec.ProteinSpecification("A", "dC", None, None), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", "dA", None, None)],
                          [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", "dA", None, None), spec.ProteinSpecification("A", "dC", None, None)]),


            OrderTestCase([spec.ProteinSpecification("A", "d", "sC", None), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", "d", "sA", None)],
                          [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", "d", "sA", None), spec.ProteinSpecification("A", "d", "sC", None)]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, "rC"), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, "rA")],
                          [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, "rA"), spec.ProteinSpecification("A", None, None, "rC")]),

            OrderTestCase([spec.ProteinSpecification("C", "dC", "sC", "rC"), spec.ProteinSpecification("B", "dB", "sB", "rB"), spec.ProteinSpecification("A", "dA", "sA", "rA")],
                          [spec.ProteinSpecification("A", "dA", "sA", "rA"), spec.ProteinSpecification("B", "dB", "sB", "rB"), spec.ProteinSpecification("C", "dC", "sC", "rC")]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, "rC"), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, None)],
                          [spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, "rC")]),

            OrderTestCase([spec.RnaSpecification("C", None, None, None), spec.RnaSpecification("B", None, None, None), spec.RnaSpecification("A", None, None, None)],
                          [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("B", None, None, None), spec.RnaSpecification("C", None, None, None)]),

            OrderTestCase([spec.RnaSpecification("A", "dC", None, None), spec.RnaSpecification("A", "dB", None, None), spec.RnaSpecification("A", "dA", None, None)],
                          [spec.RnaSpecification("A", "dA", None, None), spec.RnaSpecification("A", "dB", None, None), spec.RnaSpecification("A", "dC", None, None)]),

            OrderTestCase([spec.RnaSpecification("A", "d", "sC", None), spec.RnaSpecification("A", "d", "sB", None), spec.RnaSpecification("A", "d", "sA", None)],
                          [spec.RnaSpecification("A", "d", "sA", None), spec.RnaSpecification("A", "d", "sB", None), spec.RnaSpecification("A", "d", "sC", None)]),

            OrderTestCase([spec.RnaSpecification("A", None, None, "rC"), spec.RnaSpecification("A", None, None, "rB"), spec.RnaSpecification("A", None, None, "rA")],
                          [spec.RnaSpecification("A", None, None, "rA"), spec.RnaSpecification("A", None, None, "rB"), spec.RnaSpecification("A", None, None, "rC")]),

            OrderTestCase([spec.RnaSpecification("A", "dC", None, None), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", "dA", None, None)],
                          [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", "dA", None, None), spec.RnaSpecification("A", "dC", None, None)]),

            OrderTestCase([spec.RnaSpecification("A", "d", "sC", None), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", "d", "sA", None)],
                          [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", "d", "sA", None), spec.RnaSpecification("A", "d", "sC", None)]),

            OrderTestCase([spec.RnaSpecification("A", None, None, "rC"), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, "rA")],
                          [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, "rA"), spec.RnaSpecification("A", None, None, "rC")]),

            OrderTestCase([spec.RnaSpecification("C", "dC", "sC", "rC"), spec.RnaSpecification("B", "dB", "sB", "rB"), spec.RnaSpecification("A", "dA", "sA", "rA")],
                          [spec.RnaSpecification("A", "dA", "sA", "rA"), spec.RnaSpecification("B", "dB", "sB", "rB"), spec.RnaSpecification("C", "dC", "sC", "rC")]),

            OrderTestCase([spec.RnaSpecification("A", None, None, "rC"), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None)],
                          [spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, "rC")]),

            OrderTestCase([spec.DnaSpecification("C", None, None, None), spec.DnaSpecification("B", None, None, None), spec.DnaSpecification("A", None, None, None)],
                          [spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("B", None, None, None), spec.DnaSpecification("C", None, None, None)]),

            OrderTestCase([spec.DnaSpecification("A", "dC", None, None), spec.DnaSpecification("A", "dB", None, None), spec.DnaSpecification("A", "dA", None, None)],
                          [spec.DnaSpecification("A", "dA", None, None), spec.DnaSpecification("A", "dB", None, None), spec.DnaSpecification("A", "dC", None, None)]),

            OrderTestCase([spec.DnaSpecification("A", "d", "sC", None), spec.DnaSpecification("A", "d", "sB", None), spec.DnaSpecification("A", "d", "sA", None)],
                          [spec.DnaSpecification("A", "d", "sA", None), spec.DnaSpecification("A", "d", "sB", None), spec.DnaSpecification("A", "d", "sC", None)]),

            OrderTestCase([spec.DnaSpecification("A", None, None, "rC"), spec.DnaSpecification("A", None, None, "rB"), spec.DnaSpecification("A", None, None, "rA")],
                          [spec.DnaSpecification("A", None, None, "rA"), spec.DnaSpecification("A", None, None, "rB"), spec.DnaSpecification("A", None, None, "rC")]),

            OrderTestCase([spec.DnaSpecification("A", "dC", None, None), spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", "dA", None, None)],
                          [spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", "dA", None, None), spec.DnaSpecification("A", "dC", None, None)]),

            OrderTestCase([spec.DnaSpecification("A", "d", "sC", None), spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", "d", "sA", None)],
                          [spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", "d", "sA", None), spec.DnaSpecification("A", "d", "sC", None)]),

            OrderTestCase([spec.DnaSpecification("A", None, None, "rC"), spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", None, None, "rA")],
                          [spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", None, None, "rA"), spec.DnaSpecification("A", None, None, "rC")]),

            OrderTestCase([spec.DnaSpecification("C", "dC", "sC", "rC"), spec.DnaSpecification("B", "dB", "sB", "rB"), spec.DnaSpecification("A", "dA", "sA", "rA")],
                          [spec.DnaSpecification("A", "dA", "sA", "rA"), spec.DnaSpecification("B", "dB", "sB", "rB"), spec.DnaSpecification("C", "dC", "sC", "rC")]),

            OrderTestCase([spec.DnaSpecification("A", None, None, "rC"), spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", None, None, None)],
                          [spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", None, None, None), spec.DnaSpecification("A", None, None, "rC")])
            ]


@pytest.fixture
def the_case_class_ordering():
    return [
            OrderTestCase([spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None)],
                          [spec.RnaSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None), spec.DnaSpecification("A", None, None, None)],
                          [spec.DnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", "dA", None, None), spec.RnaSpecification("A", "dB", None, None)],
                          [spec.RnaSpecification("A", "dB", None, None), spec.ProteinSpecification("A", "dA", None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", "d", "sA", None), spec.RnaSpecification("A", "d", "sB", None)],
                          [spec.RnaSpecification("A", "d", "sB", None), spec.ProteinSpecification("A", "d", "sA", None)]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, "rA"), spec.RnaSpecification("A", None, None, "rB")],
                          [spec.RnaSpecification("A", None, None, "rB"), spec.ProteinSpecification("A", None, None, "rA")]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", "dA", None, None)],
                          [spec.RnaSpecification("A", "dA", None, None), spec.ProteinSpecification("A", None, None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", "d", "sA", None)],
                          [spec.RnaSpecification("A", "d", "sA", None), spec.ProteinSpecification("A", None, None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, "rA")],
                          [spec.RnaSpecification("A", None, None, "rA"), spec.ProteinSpecification("A", None, None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", "dA", "sA", "rA"), spec.RnaSpecification("B", "dB", "sB", "rB")],
                          [spec.RnaSpecification("B", "dB", "sB", "rB"), spec.ProteinSpecification("A", "dA", "sA", "rA")]),

            OrderTestCase([spec.RnaSpecification("B", None, None, None), spec.ProteinSpecification("B", None, None, None)],
                          [spec.RnaSpecification("B", None, None, None), spec.ProteinSpecification("B", None, None, None)]),

            OrderTestCase([spec.ProteinSpecification("A", "dA", "sA", "rA"), spec.RnaSpecification("B", "dB", "sB", "rB"), spec.RnaSpecification("A", "dB", "sB", None)],
                          [spec.RnaSpecification("A", "dB", "sB", None), spec.RnaSpecification("B", "dB", "sB", "rB"), spec.ProteinSpecification("A", "dA", "sA", "rA")]),

            OrderTestCase([spec.ProteinSpecification("A", "dA", "sA", "rA"), spec.DnaSpecification("B", "dB", "sB", "rB"), spec.RnaSpecification("A", "dB", "sB", None)],
                          [spec.DnaSpecification("B", "dB", "sB", "rB"), spec.RnaSpecification("A", "dB", "sB", None),  spec.ProteinSpecification("A", "dA", "sA", "rA")]),

            OrderTestCase([spec.ProteinSpecification("A", None, None, None), spec.DnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None)],
                          [spec.DnaSpecification("A", None, None, None), spec.RnaSpecification("A", None, None, None), spec.ProteinSpecification("A", None, None, None)])
    ]


def test_specification_internal_sorting(the_case_specifications_ordering):
    for the_case in the_case_specifications_ordering:
        assert sorted(the_case.to_sort) == the_case.expected_order


def test_specification_sorting(the_case_class_ordering):
    for the_case in the_case_class_ordering:
        assert sorted(the_case.to_sort) == the_case.expected_order

