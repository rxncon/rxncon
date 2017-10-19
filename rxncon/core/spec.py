"""Module containing the classes Spec, ProteinSpec, MRNASpec, GeneSpec, Locus, EmptyLocus, 
LocusResolution and SpecSuffix, as well as the constructor functions spec_from_str and locus_from_str.
A Spec describes a molecule-like object (be it a Protein, a Gene or an mRNA) at a certain resolution.
It consists of a component name and a Locus, the latter describing domain / residue information, and
optionally a 'struct_index', which can be used to fix the topology of complexes."""

from collections import OrderedDict
from abc import ABC
from copy import deepcopy
from typing import Optional, MutableMapping, Type, Tuple  # pylint: disable=unused-import
from enum import Enum, unique
import re

DOMAIN_SUBDOMAIN_RESIDUE_REGEX = r'^[\w:-]+\/[\w:-]+\([\w:-]+\)$'
DOMAIN_RESIDUE_REGEX = r'^[\w:-]+\([\w:-]+\)$'
DOMAIN_SUBDOMAIN_REGEX = r'^[\w:-]+\/[\w:-]+$'
RESIDUE_REGEX = r'^\([\w:-]+\)$'
DOMAIN_REGEX = r'^[\w:-]+$'


class Spec(ABC):
    """Spec is the abstract superclass for ProteinSpec, GeneSpec and MRNASpec."""

    def __init__(self, name: str, struct_index: Optional[int], locus: 'Locus') -> None:
        self.name, self.struct_index, self.locus = name, struct_index, locus
        self._validate()

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        def struct_name(spec: Spec, spec_suffix: 'SpecSuffix') -> str:
            if spec.struct_index is not None:
                return "{0}{1}@{2}".format(spec.name, spec_suffix.value, spec.struct_index)
            else:
                return "{0}{1}".format(spec.name, spec_suffix.value)

        suffix = SPEC_TO_SUFFIX[type(self)]

        if str(self.locus):
            return '{0}_[{1}]'.format(struct_name(self, suffix), str(self.locus))
        else:
            return '{0}'.format(struct_name(self, suffix))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Spec):
            return NotImplemented
        return isinstance(other, type(self)) and self.name == other.name and self.locus == other.locus and \
               self.struct_index == other.struct_index

    def __lt__(self, other: 'Spec') -> bool:
        return str(self) < str(other)

    def clone(self) -> 'Spec':
        return deepcopy(self)

    def is_subspec_of(self, other: 'Spec') -> bool:
        """Specs have a subset / superset relation, e.g. A_[(r)] is a subset of A."""
        if self == other:
            return True

        spec_pairs = zip([self.name, self.locus.domain, self.locus.subdomain, self.locus.residue],
                         [other.name, other.locus.domain, other.locus.subdomain, other.locus.residue])

        for my_property, other_property in spec_pairs:
            if my_property and other_property and my_property != other_property:
                return False
            elif not my_property and other_property:
                return False

        return True

    def is_superspec_of(self, other: 'Spec') -> bool:
        """Specs have a subset / superset relation, e.g. A is a superset of A_[(r)]."""
        if self == other:
            return True

        return other.is_subspec_of(self)

    @property
    def is_structured(self) -> bool:
        return self.struct_index is not None

    @property
    def is_component_spec(self) -> bool:
        """A ComponentSpec is a Spec without any Locus information."""
        return self.has_resolution(LocusResolution.component)

    def with_name_suffix(self, suffix: str) -> 'Spec':
        # Assure the suffix is not reserved, i.e. 'Gene', 'mRNA'
        try:
            SpecSuffix(suffix)
            raise SyntaxError('Cannot use reserved suffix {} in Spec.with_name_suffix'.format(suffix))
        except ValueError:
            return type(self)(self.name + suffix, self.struct_index, self.locus)

    def with_struct_index(self, index: int) -> 'Spec':
        """Returns a new Spec of the same type with a given struct index."""
        return type(self)(self.name, index, self.locus)

    def with_struct_from_spec(self, other: 'Spec') -> 'Spec':
        """Returns a new Spec of the same type with a struct index matching the one carried by the Spec given."""
        assert other.struct_index is not None
        assert self.name == other.name
        return type(self)(self.name, other.struct_index, self.locus)

    def with_locus(self, locus: 'Locus') -> 'Spec':
        """Returns a new Spec of the same type with a given Locus."""
        return type(self)(self.name, self.struct_index, locus.clone())

    def with_domain(self, domain: str) -> 'Spec':
        """Returns a new Spec of the same type with its domain changed into the one given."""
        return type(self)(self.name, self.struct_index, self.locus.with_domain(domain))

    def to_non_struct_spec(self) -> 'Spec':
        """Returns a new Spec of the same type without struct index."""
        return type(self)(self.name, None, self.locus)

    def to_component_spec(self) -> 'Spec':
        """Returns a new Spec of the same type without Locus."""
        return type(self)(self.name, self.struct_index, EmptyLocus())

    def to_protein_component_spec(self) -> 'ProteinSpec':
        """Returns a new ProteinSpec with the same name, without structure index or Locus."""
        return ProteinSpec(self.name, None, EmptyLocus())

    def to_gene_component_spec(self) -> 'GeneSpec':
        """Returns a new GeneSpec with the same name, without structure index or Locus."""
        return GeneSpec(self.name, None, EmptyLocus())

    def to_mrna_component_spec(self) -> 'MRNASpec':
        """Returns a new MRNASpec with the same name, without structure index or Locus."""
        return MRNASpec(self.name, None, EmptyLocus())

    @property
    def resolution(self) -> 'LocusResolution':
        return self.locus.resolution

    def has_resolution(self, resolution: 'LocusResolution') -> bool:
        return self.resolution == resolution

    def _validate(self) -> None:
        assert self.name is not None and re.match(r'\w+', self.name)


class ProteinSpec(Spec):
    def __hash__(self) -> int:
        return hash(str(self))


class MRNASpec(Spec):
    def __hash__(self) -> int:
        return hash(str(self))


class GeneSpec(Spec):
    def __hash__(self) -> int:
        return hash(str(self))


class Locus:
    """Locus contains domain, subdomain and residue information for a Spec."""

    def __init__(self, domain: Optional[str], subdomain: Optional[str], residue: Optional[str]) -> None:
        self.domain, self.subdomain, self.residue = domain, subdomain, residue
        self.validate()

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if self.domain and self.subdomain and self.residue:
            return '{0}/{1}({2})'.format(self.domain, self.subdomain, self.residue)
        elif self.domain and not self.subdomain and self.residue:
            return '{0}({1})'.format(self.domain, self.residue)
        elif self.domain and self.subdomain and not self.residue:
            return '{0}/{1}'.format(self.domain, self.subdomain)
        elif not self.domain and not self.subdomain and self.residue:
            return '({0})'.format(self.residue)
        elif self.domain and not self.subdomain and not self.residue:
            return '{0}'.format(self.domain)
        elif not self.domain and not self.subdomain and not self.residue:
            return ''
        else:
            raise AssertionError('Unable to stringify Spec')

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Locus):
            return NotImplemented
        return self.domain == other.domain and self.subdomain == other.subdomain and self.residue == other.residue

    def __lt__(self, other):
        if not isinstance(other, Locus):
            return NotImplemented
        return str(self) < str(other)

    def validate(self) -> None:
        if self.domain:
            assert re.match(r'\w+', self.domain)
        if self.subdomain:
            assert re.match(r'\w+', self.subdomain)
            assert self.domain is not None
        if self.residue:
            assert re.match(r'\w+', self.residue)

    def clone(self) -> 'Locus':
        return deepcopy(self)

    def with_domain(self, domain: str) -> 'Locus':
        locus = self.clone()
        locus.domain = domain
        locus.validate()
        return locus

    @property
    def is_empty(self) -> bool:
        return not (self.domain or self.subdomain or self.residue)

    @property
    def resolution(self) -> 'LocusResolution':
        if not self.domain and not self.subdomain and not self.residue:
            return LocusResolution.component
        elif self.domain and not self.subdomain and not self.residue:
            return LocusResolution.domain
        elif self.domain and self.subdomain and not self.residue:
            return LocusResolution.subdomain
        elif self.residue is not None:
            return LocusResolution.residue
        else:
            raise AssertionError('Inconsistent resolution for Spec {}'.format(str(self)))


def EmptyLocus() -> Locus:  # pylint: disable=invalid-name
    return Locus(None, None, None)


class LocusResolution(Enum):
    component = 'component'
    domain = 'domain'
    subdomain = 'subdomain'
    residue = 'residue'

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, LocusResolution):
            return NotImplemented

        if self == LocusResolution.component:
            return other != LocusResolution.component
        elif self == LocusResolution.domain:
            return other in [LocusResolution.subdomain, LocusResolution.residue]
        elif self == LocusResolution.subdomain:
            return other == LocusResolution.residue
        else:
            return False


@unique
class SpecSuffix(Enum):
    mrna = 'mRNA'
    dna = 'Gene'
    protein = ''


SUFFIX_TO_SPEC = OrderedDict(
    [
        (SpecSuffix.mrna, MRNASpec),
        (SpecSuffix.dna, GeneSpec),
        (SpecSuffix.protein, ProteinSpec)
    ]
)  # type: MutableMapping[SpecSuffix, Type[Spec]]

SPEC_TO_SUFFIX = OrderedDict((k, v) for v, k in SUFFIX_TO_SPEC.items())  # type: MutableMapping[Type[Spec], SpecSuffix]


def locus_from_str(locus_str: str) -> Locus:
    """Returns a Locus object parsed from the given string, raises SyntaxError if not parsable."""
    def locus_items_from_str(full_locus_str: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        domain, subdomain, residue = None, None, None  # type: Optional[str], Optional[str], Optional[str]

        if re.match(DOMAIN_SUBDOMAIN_RESIDUE_REGEX, full_locus_str):
            domain = full_locus_str.split('/')[0]
            subdomain = full_locus_str.split('/')[1].split('(')[0]
            residue = full_locus_str.split('/')[1].split('(')[1].strip(')')
        elif re.match(DOMAIN_RESIDUE_REGEX, full_locus_str):
            domain = full_locus_str.split('(')[0]
            subdomain = None
            residue = full_locus_str.split('(')[1].strip(')')
        elif re.match(DOMAIN_SUBDOMAIN_REGEX, full_locus_str):
            domain = full_locus_str.split('/')[0]
            subdomain = full_locus_str.split('/')[1]
            residue = None
        elif re.match(RESIDUE_REGEX, full_locus_str):
            domain = None
            subdomain = None
            residue = full_locus_str.strip('()')
        elif re.match(DOMAIN_REGEX, full_locus_str):
            domain = full_locus_str
            subdomain = None
            residue = None
        else:
            raise SyntaxError('Could not parse locus string {}'.format(full_locus_str))

        return domain, subdomain, residue

    return Locus(*locus_items_from_str(locus_str.strip('[]')))


def spec_from_str(spec_str: str) -> Spec:
    """Returns a Spec object parsed from the given string, raises SyntaxError if not parsable."""
    def spec_from_suffixed_name_and_locus(name: str, struct_index: Optional[int],
                                          locus: Locus) -> Spec:  # pylint: disable=invalid-name
        for suffix in SUFFIX_TO_SPEC:
            if name.endswith(suffix.value):
                name = name[:len(name) - len(suffix.value)]
                return SUFFIX_TO_SPEC[suffix](name, struct_index, locus)
            elif name.lower().endswith(suffix.value.lower()):
                raise SyntaxError(
                    'Please use correct capitalization \'{}\' in component \'{}\'.'.format(suffix.value, name))

        raise SyntaxError('Could not parse spec component_name {}'.format(name))

    if not re.match(r'[A-Za-z][A-Za-z0-9]*(@[\d]+)*(_\[[[A-Za-z0-9]/\(\)]+\])*', spec_str):
        raise SyntaxError('Spec str {} does not match validating regex.'.format(spec_str))

    struct_index = None
    items = spec_str.split('_', maxsplit=1)

    if not re.match(r'[a-zA-Z]', items[0]):
        raise SyntaxError('Spec has to start with letter character.')
    elif '@' in items[0]:
        name, struct_index_str = items[0].split('@')
        struct_index = int(struct_index_str)
    else:
        name = items[0]

    if len(items) == 1:
        return spec_from_suffixed_name_and_locus(name, struct_index, EmptyLocus())
    elif len(items) == 2:
        return spec_from_suffixed_name_and_locus(name, struct_index, locus_from_str(items[1]))
    else:
        raise SyntaxError('Could not parse spec string {}'.format(spec_str))
