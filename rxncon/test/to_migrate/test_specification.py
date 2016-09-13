import pytest
from collections import namedtuple
import rxncon.core.specification as spec

OrderTestCase = namedtuple('OrderTestCase', ['to_sort', 'expected_order'])
@pytest.fixture
def the_case_specifications_ordering():
    return [
            OrderTestCase([spec.ProteinSpec("C", None, spec.Domain(None, None, None)), spec.ProteinSpec("B", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None))],
                          [spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("B", None, spec.Domain(None, None, None)), spec.ProteinSpec("C", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("dC", None, None)), spec.ProteinSpec("A", None, spec.Domain("dB", None, None)), spec.ProteinSpec("A", None, spec.Domain("dA", None, None))],
                          [spec.ProteinSpec("A", None, spec.Domain("dA", None, None)), spec.ProteinSpec("A", None, spec.Domain("dB", None, None)), spec.ProteinSpec("A", None, spec.Domain("dC", None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("d", "sC", None)), spec.ProteinSpec("A", None, spec.Domain("d", "sB", None)), spec.ProteinSpec("A", None, spec.Domain("d", "sA", None))],
                          [spec.ProteinSpec("A", None, spec.Domain("d", "sA", None)), spec.ProteinSpec("A", None, spec.Domain("d", "sB", None)), spec.ProteinSpec("A", None, spec.Domain("d", "sC", None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, "rC")), spec.ProteinSpec("A", None, spec.Domain(None, None, "rB")), spec.ProteinSpec("A", None, spec.Domain(None, None, "rA"))],
                          [spec.ProteinSpec("A", None, spec.Domain(None, None, "rA")), spec.ProteinSpec("A", None, spec.Domain(None, None, "rB")), spec.ProteinSpec("A", None, spec.Domain(None, None, "rC"))]),


            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("dC", None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain("dA", None, None))],
                          [spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain("dA", None, None)), spec.ProteinSpec("A", None, spec.Domain("dC", None, None))]),


            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("d", "sC", None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain("d", "sA", None))],
                          [spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain("d", "sA", None)), spec.ProteinSpec("A", None, spec.Domain("d", "sC", None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, "rC")), spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, "rA"))],
                          [spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, "rA")), spec.ProteinSpec("A", None, spec.Domain(None, None, "rC"))]),

            OrderTestCase([spec.ProteinSpec("C", None, spec.Domain("dC", "sC", "rC")), spec.ProteinSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.ProteinSpec("A", None, spec.Domain("dA", "sA", "rA"))],
                          [spec.ProteinSpec("A", None, spec.Domain("dA", "sA", "rA")), spec.ProteinSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.ProteinSpec("C", None, spec.Domain("dC", "sC", "rC"))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, "rC")), spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None))],
                          [spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpec("C", None, spec.Domain(None, None, None)), spec.RnaSpec("B", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("B", None, spec.Domain(None, None, None)), spec.RnaSpec("C", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.RnaSpec("A", None, spec.Domain("dC", None, None)), spec.RnaSpec("A", None, spec.Domain("dB", None, None)), spec.RnaSpec("A", None, spec.Domain("dA", None, None))],
                          [spec.RnaSpec("A", None, spec.Domain("dA", None, None)), spec.RnaSpec("A", None, spec.Domain("dB", None, None)), spec.RnaSpec("A", None, spec.Domain("dC", None, None))]),

            OrderTestCase([spec.RnaSpec("A", None, spec.Domain("d", "sC", None)), spec.RnaSpec("A", None, spec.Domain("d", "sB", None)), spec.RnaSpec("A", None, spec.Domain("d", "sA", None))],
                          [spec.RnaSpec("A", None, spec.Domain("d", "sA", None)), spec.RnaSpec("A", None, spec.Domain("d", "sB", None)), spec.RnaSpec("A", None, spec.Domain("d", "sC", None))]),

            OrderTestCase([spec.RnaSpec("A", None, spec.Domain(None, None, "rC")), spec.RnaSpec("A", None, spec.Domain(None, None, "rB")), spec.RnaSpec("A", None, spec.Domain(None, None, "rA"))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, "rA")), spec.RnaSpec("A", None, spec.Domain(None, None, "rB")), spec.RnaSpec("A", None, spec.Domain(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpec("A", None, spec.Domain("dC", None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain("dA", None, None))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain("dA", None, None)), spec.RnaSpec("A", None, spec.Domain("dC", None, None))]),

            OrderTestCase([spec.RnaSpec("A", None, spec.Domain("d", "sC", None)), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain("d", "sA", None))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain("d", "sA", None)), spec.RnaSpec("A", None, spec.Domain("d", "sC", None))]),

            OrderTestCase([spec.RnaSpec("A", None, spec.Domain(None, None, "rC")), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, "rA"))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, "rA")), spec.RnaSpec("A", None, spec.Domain(None, None, "rC"))]),

            OrderTestCase([spec.RnaSpec("C", None, spec.Domain("dC", "sC", "rC")), spec.RnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.RnaSpec("A", None, spec.Domain("dA", "sA", "rA"))],
                          [spec.RnaSpec("A", None, spec.Domain("dA", "sA", "rA")), spec.RnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.RnaSpec("C", None, spec.Domain("dC", "sC", "rC"))]),

            OrderTestCase([spec.RnaSpec("A", None, spec.Domain(None, None, "rC")), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpec("C", None, spec.Domain(None, None, None)), spec.DnaSpec("B", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, None))],
                          [spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("B", None, spec.Domain(None, None, None)), spec.DnaSpec("C", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.DnaSpec("A", None, spec.Domain("dC", None, None)), spec.DnaSpec("A", None, spec.Domain("dB", None, None)), spec.DnaSpec("A", None, spec.Domain("dA", None, None))],
                          [spec.DnaSpec("A", None, spec.Domain("dA", None, None)), spec.DnaSpec("A", None, spec.Domain("dB", None, None)), spec.DnaSpec("A", None, spec.Domain("dC", None, None))]),

            OrderTestCase([spec.DnaSpec("A", None, spec.Domain("d", "sC", None)), spec.DnaSpec("A", None, spec.Domain("d", "sB", None)), spec.DnaSpec("A", None, spec.Domain("d", "sA", None))],
                          [spec.DnaSpec("A", None, spec.Domain("d", "sA", None)), spec.DnaSpec("A", None, spec.Domain("d", "sB", None)), spec.DnaSpec("A", None, spec.Domain("d", "sC", None))]),

            OrderTestCase([spec.DnaSpec("A", None, spec.Domain(None, None, "rC")), spec.DnaSpec("A", None, spec.Domain(None, None, "rB")), spec.DnaSpec("A", None, spec.Domain(None, None, "rA"))],
                          [spec.DnaSpec("A", None, spec.Domain(None, None, "rA")), spec.DnaSpec("A", None, spec.Domain(None, None, "rB")), spec.DnaSpec("A", None, spec.Domain(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpec("A", None, spec.Domain("dC", None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain("dA", None, None))],
                          [spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain("dA", None, None)), spec.DnaSpec("A", None, spec.Domain("dC", None, None))]),

            OrderTestCase([spec.DnaSpec("A", None, spec.Domain("d", "sC", None)), spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain("d", "sA", None))],
                          [spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain("d", "sA", None)), spec.DnaSpec("A", None, spec.Domain("d", "sC", None))]),

            OrderTestCase([spec.DnaSpec("A", None, spec.Domain(None, None, "rC")), spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, "rA"))],
                          [spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, "rA")), spec.DnaSpec("A", None, spec.Domain(None, None, "rC"))]),

            OrderTestCase([spec.DnaSpec("C", None, spec.Domain("dC", "sC", "rC")), spec.DnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.DnaSpec("A", None, spec.Domain("dA", "sA", "rA"))],
                          [spec.DnaSpec("A", None, spec.Domain("dA", "sA", "rA")), spec.DnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.DnaSpec("C", None, spec.Domain("dC", "sC", "rC"))]),

            OrderTestCase([spec.DnaSpec("A", None, spec.Domain(None, None, "rC")), spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, None))],
                          [spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, "rC"))])
            ]


@pytest.fixture
def the_case_class_ordering():
    return [
            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, None))],
                          [spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("dA", None, None)), spec.RnaSpec("A", None, spec.Domain("dB", None, None))],
                          [spec.RnaSpec("A", None, spec.Domain("dB", None, None)), spec.ProteinSpec("A", None, spec.Domain("dA", None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("d", "sA", None)), spec.RnaSpec("A", None, spec.Domain("d", "sB", None))],
                          [spec.RnaSpec("A", None, spec.Domain("d", "sB", None)), spec.ProteinSpec("A", None, spec.Domain("d", "sA", None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, "rA")), spec.RnaSpec("A", None, spec.Domain(None, None, "rB"))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, "rB")), spec.ProteinSpec("A", None, spec.Domain(None, None, "rA"))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain("dA", None, None))],
                          [spec.RnaSpec("A", None, spec.Domain("dA", None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain("d", "sA", None))],
                          [spec.RnaSpec("A", None, spec.Domain("d", "sA", None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, "rA"))],
                          [spec.RnaSpec("A", None, spec.Domain(None, None, "rA")), spec.ProteinSpec("A", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("dA", "sA", "rA")), spec.RnaSpec("B", None, spec.Domain("dB", "sB", "rB"))],
                          [spec.RnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.ProteinSpec("A", None, spec.Domain("dA", "sA", "rA"))]),

            OrderTestCase([spec.RnaSpec("B", None, spec.Domain(None, None, None)), spec.ProteinSpec("B", None, spec.Domain(None, None, None))],
                          [spec.RnaSpec("B", None, spec.Domain(None, None, None)), spec.ProteinSpec("B", None, spec.Domain(None, None, None))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("dA", "sA", "rA")), spec.RnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.RnaSpec("A", None, spec.Domain("dB", "sB", None))],
                          [spec.RnaSpec("A", None, spec.Domain("dB", "sB", None)), spec.RnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.ProteinSpec("A", None, spec.Domain("dA", "sA", "rA"))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain("dA", "sA", "rA")), spec.DnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.RnaSpec("A", None, spec.Domain("dB", "sB", None))],
                          [spec.DnaSpec("B", None, spec.Domain("dB", "sB", "rB")), spec.RnaSpec("A", None, spec.Domain("dB", "sB", None)), spec.ProteinSpec("A", None, spec.Domain("dA", "sA", "rA"))]),

            OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None))],
                          [spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None))]),

        OrderTestCase([spec.ProteinSpec("A", None, spec.Domain(None, None, None)), spec.EmptySpec("0", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None))],
                      [spec.EmptySpec("0", None, spec.Domain(None, None, None)), spec.DnaSpec("A", None, spec.Domain(None, None, None)), spec.RnaSpec("A", None, spec.Domain(None, None, None)), spec.ProteinSpec("A", None, spec.Domain(None, None, None))])
    ]

StringTestCase = namedtuple('StringTestCase', ['specification', 'expected_string'])

@pytest.fixture
def the_case_string_generation():
    return [StringTestCase(spec.ProteinSpec("C", None, spec.Domain(None, None, None)),
                           'C'),
            StringTestCase(spec.ProteinSpec("C", None, spec.Domain('d', 's', None)),
                          'C_[d/s]'),
            StringTestCase(spec.ProteinSpec("C", None, spec.Domain(None, None, 'r')),
                          'C_[(r)]'),
            StringTestCase(spec.ProteinSpec("C", None, spec.Domain('d', None, 'r')),
                           'C_[d(r)]'),
            StringTestCase(spec.ProteinSpec("C", None, spec.Domain('d', 's', 'r')),
                           'C_[d/s(r)]'),
            StringTestCase(spec.ProteinSpec("C", 0, spec.Domain('d', 's', 'r')),
                           'C@0_[d/s(r)]'),

            StringTestCase(spec.EmptySpec("0", None, spec.Domain(None, None, None)),
                           '0'),

            StringTestCase(spec.DnaSpec("C", None, spec.Domain(None, None, None)),
                           'Cgene'),
            StringTestCase(spec.DnaSpec("C", None, spec.Domain('d', 's', None)),
                           'Cgene_[d/s]'),
            StringTestCase(spec.DnaSpec("C", None, spec.Domain(None, None, 'r')),
                           'Cgene_[(r)]'),
            StringTestCase(spec.DnaSpec("C", None, spec.Domain('d', None, 'r')),
                           'Cgene_[d(r)]'),
            StringTestCase(spec.DnaSpec("C", None, spec.Domain('d', 's', 'r')),
                           'Cgene_[d/s(r)]'),
            StringTestCase(spec.DnaSpec("C", 0, spec.Domain('d', 's', 'r')),
                           'Cgene@0_[d/s(r)]'),

            StringTestCase(spec.RnaSpec("C", None, spec.Domain(None, None, None)),
                           'CmRNA'),
            StringTestCase(spec.RnaSpec("C", None, spec.Domain('d', 's', None)),
                           'CmRNA_[d/s]'),
            StringTestCase(spec.RnaSpec("C", None, spec.Domain(None, None, 'r')),
                           'CmRNA_[(r)]'),
            StringTestCase(spec.RnaSpec("C", None, spec.Domain('d', None, 'r')),
                           'CmRNA_[d(r)]'),
            StringTestCase(spec.RnaSpec("C", None, spec.Domain('d', 's', 'r')),
                           'CmRNA_[d/s(r)]'),
            StringTestCase(spec.RnaSpec("C", 0, spec.Domain('d', 's', 'r')),
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


