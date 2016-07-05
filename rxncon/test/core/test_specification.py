import pytest
from collections import namedtuple
import rxncon.core.specification as spec

OrderTestCase = namedtuple('OrderTestCase', ['to_sort', 'expected_order'])
@pytest.fixture
def the_case_specifications_ordering():
    return [
            OrderTestCase([spec.ProteinSpecification("C", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("B", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("B", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("C", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("dC", None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dB", None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", None, None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dB", None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dC", None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sC", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sB", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sA", None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sA", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sB", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sC", None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rB")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rA"))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rA")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rB")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))]),


            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("dC", None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", None, None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dC", None, None))]),


            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sC", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sA", None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sA", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sC", None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rA"))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rA")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))]),

            OrderTestCase([spec.ProteinSpecification("C", 0, spec.DomainDefinition("dC", "sC", "rC")), spec.ProteinSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA"))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA")), spec.ProteinSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.ProteinSpecification("C", 0, spec.DomainDefinition("dC", "sC", "rC"))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("C", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("B", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("B", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("C", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainDefinition("dC", None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dB", None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dA", None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition("dA", None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dB", None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dC", None, None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sC", None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sB", None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sB", None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sC", None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rB")), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA"))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA")), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rB")), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainDefinition("dC", None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dA", None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dA", None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dC", None, None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sC", None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sC", None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA"))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA")), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("C", 0, spec.DomainDefinition("dC", "sC", "rC")), spec.RnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.RnaSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA"))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA")), spec.RnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.RnaSpecification("C", 0, spec.DomainDefinition("dC", "sC", "rC"))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("C", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("B", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("B", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("C", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainDefinition("dC", None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("dB", None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("dA", None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition("dA", None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("dB", None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("dC", None, None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sC", None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sB", None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sB", None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sC", None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rB")), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA"))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA")), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rB")), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainDefinition("dC", None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("dA", None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("dA", None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("dC", None, None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sC", None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None)), spec.DnaSpecification("A", 0, spec.DomainDefinition("d", "sC", None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA"))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA")), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("C", 0, spec.DomainDefinition("dC", "sC", "rC")), spec.DnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.DnaSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA"))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA")), spec.DnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.DnaSpecification("C", 0, spec.DomainDefinition("dC", "sC", "rC"))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC")), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, "rC"))])
            ]


@pytest.fixture
def the_case_class_ordering():
    return [
            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dB", None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition("dB", None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sA", None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sB", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sB", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("d", "sA", None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rA")), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rB"))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rB")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("dA", None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition("dA", None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition("d", "sA", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA"))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, "rA")), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA")), spec.RnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB"))],
                          [spec.RnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA"))]),

            OrderTestCase([spec.RnaSpecification("B", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("B", 0, spec.DomainDefinition(None, None, None))],
                          [spec.RnaSpecification("B", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("B", 0, spec.DomainDefinition(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA")), spec.RnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.RnaSpecification("A", 0, spec.DomainDefinition("dB", "sB", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainDefinition("dB", "sB", None)), spec.RnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA")), spec.DnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.RnaSpecification("A", 0, spec.DomainDefinition("dB", "sB", None))],
                          [spec.DnaSpecification("B", 0, spec.DomainDefinition("dB", "sB", "rB")), spec.RnaSpecification("A", 0, spec.DomainDefinition("dB", "sB", None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition("dA", "sA", "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainDefinition(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainDefinition(None, None, None))])
    ]

StringTestCase = namedtuple('StringTestCase', ['specification', 'expected_string'])

@pytest.fixture
def the_case_string_generation():
    return [StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainDefinition(None, None, None)),
                           'C@0'),
            StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainDefinition('d', 's', None)),
                          'C@0_[d/s]'),
            StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainDefinition(None, None, 'r')),
                          'C@0_[(r)]'),
            StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainDefinition('d', None, 'r')),
                           'C@0_[d(r)]'),
            StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainDefinition('d', 's', 'r')),
                           'C@0_[d/s(r)]'),

            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainDefinition(None, None, None)),
                           'Cgene@0'),
            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainDefinition('d', 's', None)),
                           'Cgene@0_[d/s]'),
            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainDefinition(None, None, 'r')),
                           'Cgene@0_[(r)]'),
            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainDefinition('d', None, 'r')),
                           'Cgene@0_[d(r)]'),
            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainDefinition('d', 's', 'r')),
                           'Cgene@0_[d/s(r)]'),

            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainDefinition(None, None, None)),
                           'CmRNA@0'),
            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainDefinition('d', 's', None)),
                           'CmRNA@0_[d/s]'),
            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainDefinition(None, None, 'r')),
                           'CmRNA@0_[(r)]'),
            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainDefinition('d', None, 'r')),
                           'CmRNA@0_[d(r)]'),
            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainDefinition('d', 's', 'r')),
                           'CmRNA@0_[d/s(r)]'),

            ]


def test_specification_internal_sorting(the_case_specifications_ordering):
    for the_case in the_case_specifications_ordering:
        assert sorted(the_case.to_sort) == the_case.expected_order


def test_specification_sorting(the_case_class_ordering):
    for the_case in the_case_class_ordering:
        assert sorted(the_case.to_sort) == the_case.expected_order


def test_string_generation(the_case_string_generation):
    for the_case in the_case_string_generation:
        assert str(the_case.specification) == the_case.expected_string


