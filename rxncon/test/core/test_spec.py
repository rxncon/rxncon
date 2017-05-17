import pytest
from collections import namedtuple

from rxncon.core.spec import Locus, GeneSpec, MRNASpec, ProteinSpec, locus_from_str, LocusResolution, spec_from_str


def test_loci() -> None:
    LocusTestCase = namedtuple('LocusTestCase',
                               ['locus_str', 'expected_domain', 'expected_subdomain', 'expected_residue'])

    locus_cases = [
        LocusTestCase('d', 'd', None, None),
        LocusTestCase('d/s', 'd', 's', None),
        LocusTestCase('d/s(r)', 'd', 's', 'r'),
        LocusTestCase('(r)', None, None, 'r')
    ]

    for locus_case in locus_cases:
        assert locus_from_str(locus_case.locus_str) == \
               Locus(locus_case.expected_domain, locus_case.expected_subdomain, locus_case.expected_residue)


def test_unstructured_specs() -> None:
    protein_spec = spec_from_str('A_[dd/ss(rr)')
    # Protein
    assert isinstance(protein_spec, ProteinSpec)
    assert protein_spec.has_resolution(LocusResolution.residue)
    assert protein_spec.to_component_spec().has_resolution(LocusResolution.component)
    # DNA
    assert isinstance(protein_spec.to_gene_component_spec(), GeneSpec)
    assert protein_spec.to_gene_component_spec().has_resolution(LocusResolution.component)
    # mRNA
    assert isinstance(protein_spec.to_mrna_component_spec(), MRNASpec)
    assert protein_spec.to_mrna_component_spec().has_resolution(LocusResolution.component)

    with pytest.raises(SyntaxError):
        spec_from_str('0')


def test_structured_specs() -> None:
    protein_spec = spec_from_str('A@0_[dd/ss(rr)')
    # Protein
    assert isinstance(protein_spec, ProteinSpec)
    assert protein_spec.has_resolution(LocusResolution.residue)
    assert protein_spec.to_component_spec().has_resolution(LocusResolution.component)
    assert protein_spec.struct_index == 0
    # DNA
    assert isinstance(protein_spec.to_gene_component_spec(), GeneSpec)
    assert protein_spec.to_gene_component_spec().has_resolution(LocusResolution.component)
    assert not protein_spec.to_gene_component_spec().struct_index
    # mRNA
    assert isinstance(protein_spec.to_mrna_component_spec(), MRNASpec)
    assert protein_spec.to_mrna_component_spec().has_resolution(LocusResolution.component)
    assert not protein_spec.to_mrna_component_spec().struct_index

    with pytest.raises(SyntaxError):
        spec_from_str('0@1')


def test_super_sub() -> None:
    assert spec_from_str('A_[x]').is_subspec_of(spec_from_str('A'))
    assert not spec_from_str('A').is_subspec_of(spec_from_str('B'))

    assert spec_from_str('A_[dom]').is_subspec_of(spec_from_str('A'))
    assert not spec_from_str('A').is_subspec_of(spec_from_str('A_[dom]'))


def test_equality() -> None:
    non_struct = spec_from_str('A_[d/s(r)]')
    struct     = spec_from_str('A@5_[d/s(r)]')

    assert struct.to_non_struct_spec() == non_struct
    assert not struct == non_struct
    assert not non_struct == struct


def test_colons_dashes() -> None:
    assert spec_from_str('A_[d:bla]').locus.domain == 'd:bla'
    assert spec_from_str('A_[d-bla]').locus.domain == 'd-bla'
    assert spec_from_str('A_[(r:bla)]').locus.residue == 'r:bla'
    assert spec_from_str('A_[(r-bla)]').locus.residue == 'r-bla'
