import pytest
from collections import namedtuple
import rxncon.core.specification as spec

OrderTestCase = namedtuple('OrderTestCase', ['to_sort', 'expected_order'])
@pytest.fixture
def the_case_specifications_ordering():
    return [
            OrderTestCase([spec.ProteinSpecification("C", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("B", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.ProteinSpecification("A", 0,  spec.DomainResolution(None, None, None)), spec.ProteinSpecification("B", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("C", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("dC", None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("dB", None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", None, None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("dB", None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sC", None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sB", None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sA", None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sA", None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sB", None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rB")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rA"))],
                          [spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rA")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rB")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rC"))]),


            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("dC", None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", None, None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("dC", None, None))]),


            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sC", None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sA", None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sA", None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rA"))],
                          [spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rA")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.ProteinSpecification("C", 0, spec.DomainResolution("dC", "sC", "rC")), spec.ProteinSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA"))],
                          [spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA")), spec.ProteinSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.ProteinSpecification("C", 0, spec.DomainResolution("dC", "sC", "rC"))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("C", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("B", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("B", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("C", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainResolution("dC", None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dB", None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dA", None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution("dA", None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dB", None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sC", None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sB", None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sA", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sA", None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sB", None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rB")), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rA"))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rA")), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rB")), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainResolution("dC", None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dA", None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dA", None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sC", None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sA", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sA", None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rA"))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rA")), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpecification("C", 0, spec.DomainResolution("dC", "sC", "rC")), spec.RnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA"))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA")), spec.RnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("C", 0, spec.DomainResolution("dC", "sC", "rC"))]),

            OrderTestCase([spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("C", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("B", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("B", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("C", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainResolution("dC", None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("dB", None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("dA", None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution("dA", None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("dB", None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sC", None)), spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sB", None)), spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sA", None))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sA", None)), spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sB", None)), spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rB")), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rA"))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rA")), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rB")), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainResolution("dC", None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("dA", None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("dA", None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("dC", None, None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sC", None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sA", None))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sA", None)), spec.DnaSpecification("A", 0, spec.DomainResolution("d", "sC", None))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rA"))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rA")), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpecification("C", 0, spec.DomainResolution("dC", "sC", "rC")), spec.DnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.DnaSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA"))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA")), spec.DnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.DnaSpecification("C", 0, spec.DomainResolution("dC", "sC", "rC"))]),

            OrderTestCase([spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rC")), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, "rC"))])
            ]


@pytest.fixture
def the_case_class_ordering():
    return [
            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dB", None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution("dB", None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sA", None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sB", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sB", None)), spec.ProteinSpecification("A", 0, spec.DomainResolution("d", "sA", None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rA")), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rB"))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rB")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("dA", None, None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution("dA", None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sA", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution("d", "sA", None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rA"))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, "rA")), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA")), spec.RnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB"))],
                          [spec.RnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA"))]),

            OrderTestCase([spec.RnaSpecification("B", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("B", 0, spec.DomainResolution(None, None, None))],
                          [spec.RnaSpecification("B", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("B", 0, spec.DomainResolution(None, None, None))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA")), spec.RnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("A", 0, spec.DomainResolution("dB", "sB", None))],
                          [spec.RnaSpecification("A", 0, spec.DomainResolution("dB", "sB", None)), spec.RnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA")), spec.DnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("A", 0, spec.DomainResolution("dB", "sB", None))],
                          [spec.DnaSpecification("B", 0, spec.DomainResolution("dB", "sB", "rB")), spec.RnaSpecification("A", 0, spec.DomainResolution("dB", "sB", None)),  spec.ProteinSpecification("A", 0, spec.DomainResolution("dA", "sA", "rA"))]),

            OrderTestCase([spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None))],
                          [spec.DnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.RnaSpecification("A", 0, spec.DomainResolution(None, None, None)), spec.ProteinSpecification("A", 0, spec.DomainResolution(None, None, None))])
    ]

StringTestCase = namedtuple('StringTestCase', ['specification', 'expected_string'])

@pytest.fixture
def the_case_string_generation():
    return [StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainResolution(None, None, None)),
                           'ProteinSpecification: C@0'),
            StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainResolution('d', 's', None)),
                          'ProteinSpecification: C@0_[d/s]'),
            StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainResolution(None, None, 'r')),
                          'ProteinSpecification: C@0_[(r)]'),
            StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainResolution('d', None, 'r')),
                           'ProteinSpecification: C@0_[d(r)]'),
            StringTestCase(spec.ProteinSpecification("C", 0, spec.DomainResolution('d', 's', 'r')),
                           'ProteinSpecification: C@0_[d/s(r)]'),

            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainResolution(None, None, None)),
                           'DnaSpecification: C@0gene'),
            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainResolution('d', 's', None)),
                           'DnaSpecification: C@0gene_[d/s]'),
            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainResolution(None, None, 'r')),
                           'DnaSpecification: C@0gene_[(r)]'),
            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainResolution('d', None, 'r')),
                           'DnaSpecification: C@0gene_[d(r)]'),
            StringTestCase(spec.DnaSpecification("C", 0, spec.DomainResolution('d', 's', 'r')),
                           'DnaSpecification: C@0gene_[d/s(r)]'),

            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainResolution(None, None, None)),
                           'RnaSpecification: C@0mRNA'),
            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainResolution('d', 's', None)),
                           'RnaSpecification: C@0mRNA_[d/s]'),
            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainResolution(None, None, 'r')),
                           'RnaSpecification: C@0mRNA_[(r)]'),
            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainResolution('d', None, 'r')),
                           'RnaSpecification: C@0mRNA_[d(r)]'),
            StringTestCase(spec.RnaSpecification("C", 0, spec.DomainResolution('d', 's', 'r')),
                           'RnaSpecification: C@0mRNA_[d/s(r)]'),

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


