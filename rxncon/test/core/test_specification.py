import pytest
from collections import namedtuple
import rxncon.core.specification as spec

OrderTestCase = namedtuple('OrderTestCase', ['to_sort', 'expected_order'])
@pytest.fixture
def the_case_specifications_ordering():
    return [
            OrderTestCase([spec.ProteinSpecification("C", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("B", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("B", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("C", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("dC", None, None)), spec.ProteinSpecification("A", spec.DomainResolution("dB", None, None)), spec.ProteinSpecification("A", spec.DomainResolution("dA", None, None))],
                          [spec.ProteinSpecification("A", spec.DomainResolution("dA", None, None)), spec.ProteinSpecification("A", spec.DomainResolution("dB", None, None)), spec.ProteinSpecification("A", spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("d", "sC", None)), spec.ProteinSpecification("A", spec.DomainResolution("d", "sB", None)), spec.ProteinSpecification("A", spec.DomainResolution("d", "sA", None))],
                          [spec.ProteinSpecification("A", spec.DomainResolution("d", "sA", None)), spec.ProteinSpecification("A", spec.DomainResolution("d", "sB", None)), spec.ProteinSpecification("A", spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rC")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rB")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rA"))],
                          [spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rA")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rB")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rC"))]),


            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("dC", None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution("dA", None, None))],
                          [spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution("dA", None, None)), spec.ProteinSpecification("A", spec.DomainResolution("dC", None, None))]),


            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("d", "sC", None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution("d", "sA", None))],
                          [spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution("d", "sA", None)), spec.ProteinSpecification("A", spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rC")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rA"))],
                          [spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rA")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.ProteinSpecification("C", spec.DomainResolution("dC", "sC", "rC")), spec.ProteinSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.ProteinSpecification("A", spec.DomainResolution("dA", "sA", "rA"))],
                          [spec.ProteinSpecification("A", spec.DomainResolution("dA", "sA", "rA")), spec.ProteinSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.ProteinSpecification("C", spec.DomainResolution("dC", "sC", "rC"))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rC")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("C", spec.DomainResolution(None, None, None)), spec.RnaSpecification("B", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("B", spec.DomainResolution(None, None, None)), spec.RnaSpecification("C", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.RnaSpecification("A", spec.DomainResolution("dC", None, None)), spec.RnaSpecification("A", spec.DomainResolution("dB", None, None)), spec.RnaSpecification("A", spec.DomainResolution("dA", None, None))],
                          [spec.RnaSpecification("A", spec.DomainResolution("dA", None, None)), spec.RnaSpecification("A", spec.DomainResolution("dB", None, None)), spec.RnaSpecification("A", spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.RnaSpecification("A", spec.DomainResolution("d", "sC", None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sB", None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sA", None))],
                          [spec.RnaSpecification("A", spec.DomainResolution("d", "sA", None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sB", None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.RnaSpecification("A", spec.DomainResolution(None, None, "rC")), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rB")), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rA"))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, "rA")), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rB")), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("A", spec.DomainResolution("dC", None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution("dA", None, None))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution("dA", None, None)), spec.RnaSpecification("A", spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.RnaSpecification("A", spec.DomainResolution("d", "sC", None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sA", None))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sA", None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.RnaSpecification("A", spec.DomainResolution(None, None, "rC")), spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rA"))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rA")), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("C", spec.DomainResolution("dC", "sC", "rC")), spec.RnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("A", spec.DomainResolution("dA", "sA", "rA"))],
                          [spec.RnaSpecification("A", spec.DomainResolution("dA", "sA", "rA")), spec.RnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("C", spec.DomainResolution("dC", "sC", "rC"))]),

            OrderTestCase([spec.RnaSpecification("A", spec.DomainResolution(None, None, "rC")), spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("C", spec.DomainResolution(None, None, None)), spec.DnaSpecification("B", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("B", spec.DomainResolution(None, None, None)), spec.DnaSpecification("C", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.DnaSpecification("A", spec.DomainResolution("dC", None, None)), spec.DnaSpecification("A", spec.DomainResolution("dB", None, None)), spec.DnaSpecification("A", spec.DomainResolution("dA", None, None))],
                          [spec.DnaSpecification("A", spec.DomainResolution("dA", None, None)), spec.DnaSpecification("A", spec.DomainResolution("dB", None, None)), spec.DnaSpecification("A", spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.DnaSpecification("A", spec.DomainResolution("d", "sC", None)), spec.DnaSpecification("A", spec.DomainResolution("d", "sB", None)), spec.DnaSpecification("A", spec.DomainResolution("d", "sA", None))],
                          [spec.DnaSpecification("A", spec.DomainResolution("d", "sA", None)), spec.DnaSpecification("A", spec.DomainResolution("d", "sB", None)), spec.DnaSpecification("A", spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.DnaSpecification("A", spec.DomainResolution(None, None, "rC")), spec.DnaSpecification("A", spec.DomainResolution(None, None, "rB")), spec.DnaSpecification("A", spec.DomainResolution(None, None, "rA"))],
                          [spec.DnaSpecification("A", spec.DomainResolution(None, None, "rA")), spec.DnaSpecification("A", spec.DomainResolution(None, None, "rB")), spec.DnaSpecification("A", spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("A", spec.DomainResolution("dC", None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution("dA", None, None))],
                          [spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution("dA", None, None)), spec.DnaSpecification("A", spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.DnaSpecification("A", spec.DomainResolution("d", "sC", None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution("d", "sA", None))],
                          [spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution("d", "sA", None)), spec.DnaSpecification("A", spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.DnaSpecification("A", spec.DomainResolution(None, None, "rC")), spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, "rA"))],
                          [spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, "rA")), spec.DnaSpecification("A", spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("C", spec.DomainResolution("dC", "sC", "rC")), spec.DnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.DnaSpecification("A", spec.DomainResolution("dA", "sA", "rA"))],
                          [spec.DnaSpecification("A", spec.DomainResolution("dA", "sA", "rA")), spec.DnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.DnaSpecification("C", spec.DomainResolution("dC", "sC", "rC"))]),

            OrderTestCase([spec.DnaSpecification("A", spec.DomainResolution(None, None, "rC")), spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, "rC"))])
            ]


@pytest.fixture
def the_case_class_ordering():
    return [
            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("dA", None, None)), spec.RnaSpecification("A", spec.DomainResolution("dB", None, None))],
                          [spec.RnaSpecification("A", spec.DomainResolution("dB", None, None)), spec.ProteinSpecification("A", spec.DomainResolution("dA", None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("d", "sA", None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sB", None))],
                          [spec.RnaSpecification("A", spec.DomainResolution("d", "sB", None)), spec.ProteinSpecification("A", spec.DomainResolution("d", "sA", None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rA")), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rB"))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, "rB")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution("dA", None, None))],
                          [spec.RnaSpecification("A", spec.DomainResolution("dA", None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution("d", "sA", None))],
                          [spec.RnaSpecification("A", spec.DomainResolution("d", "sA", None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, "rA"))],
                          [spec.RnaSpecification("A", spec.DomainResolution(None, None, "rA")), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("dA", "sA", "rA")), spec.RnaSpecification("B", spec.DomainResolution("dB", "sB", "rB"))],
                          [spec.RnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.ProteinSpecification("A", spec.DomainResolution("dA", "sA", "rA"))]),

            OrderTestCase([spec.RnaSpecification("B", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("B", spec.DomainResolution(None, None, None))],
                          [spec.RnaSpecification("B", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("B", spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("dA", "sA", "rA")), spec.RnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("A", spec.DomainResolution("dB", "sB", None))],
                          [spec.RnaSpecification("A", spec.DomainResolution("dB", "sB", None)), spec.RnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.ProteinSpecification("A", spec.DomainResolution("dA", "sA", "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution("dA", "sA", "rA")), spec.DnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("A", spec.DomainResolution("dB", "sB", None))],
                          [spec.DnaSpecification("B", spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("A", spec.DomainResolution("dB", "sB", None)),  spec.ProteinSpecification("A", spec.DomainResolution("dA", "sA", "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None))],
                          [spec.DnaSpecification("A", spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", spec.DomainResolution(None, None, None))])
    ]


def test_specification_internal_sorting(the_case_specifications_ordering):
    for the_case in the_case_specifications_ordering:
        assert sorted(the_case.to_sort) == the_case.expected_order


def test_specification_sorting(the_case_class_ordering):
    for the_case in the_case_class_ordering:
        assert sorted(the_case.to_sort) == the_case.expected_order

