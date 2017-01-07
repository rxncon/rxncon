from collections import OrderedDict
from abc import ABC
from typing import Optional
from enum import Enum, unique
import re

from rxncon.util.utils import OrderedEnum


class Spec(ABC):
    def __init__(self, component_name: str, struct_index: Optional[int], locus: 'Locus'):
        self.component_name, self.struct_index, self.locus = component_name, struct_index, locus
        self._validate()

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        def struct_name(spec: Spec, suffix: 'SpecSuffix'):
            if spec.struct_index is not None:
                return "{0}{1}@{2}".format(spec.component_name, suffix.value, spec.struct_index)
            else:
                return "{0}{1}".format(spec.component_name, suffix.value)

        suffix = spec_to_suffix[type(self)]

        if str(self.locus):
            return '{0}_[{1}]'.format(struct_name(self, suffix), str(self.locus))
        else:
            return '{0}'.format(struct_name(self, suffix))

    def __eq__(self, other: 'Spec') -> bool:
        return isinstance(other, type(self)) and self.component_name == other.component_name and self.locus == other.locus \
            and self.struct_index == other.struct_index

    def __lt__(self, other: 'Spec') -> bool:
        return str(self) < str(other)

    def is_subspec_of(self, other: 'Spec') -> bool:
        if self == other:
            return True

        spec_pairs = zip([self.component_name, self.locus.domain, self.locus.subdomain, self.locus.residue],
                         [other.component_name, other.locus.domain, other.locus.subdomain, other.locus.residue])

        for my_property, other_property in spec_pairs:
            if my_property and other_property and my_property != other_property:
                return False
            elif not my_property and other_property:
                return False

        return True

    def is_superspec_of(self, other: 'Spec') -> bool:
        if self == other:
            return True

        return other.is_subspec_of(self)

    @property
    def is_structured(self):
        return self.struct_index is not None

    @property
    def is_component_spec(self) -> bool:
        return self.has_resolution(LocusResolution.component)

    def with_struct_index(self, index: int) -> 'Spec':
        return type(self)(self.component_name, index, self.locus)

    def with_struct_from_spec(self, other: 'Spec') -> 'Spec':
        assert other.struct_index is not None
        assert self.component_name == other.component_name
        return type(self)(self.component_name, other.struct_index, self.locus)

    def to_non_struct_spec(self):
        return type(self)(self.component_name, None, self.locus)

    def to_component_spec(self) -> 'Spec':
        return type(self)(self.component_name, self.struct_index, EmptyLocus())

    def to_protein_component_spec(self) -> 'ProteinSpec':
        return ProteinSpec(self.component_name, None, EmptyLocus())

    def to_dna_component_spec(self) -> 'GeneSpec':
        return GeneSpec(self.component_name, None, EmptyLocus())

    def to_mrna_component_spec(self) -> 'MRNASpec':
        return MRNASpec(self.component_name, None, EmptyLocus())

    @property
    def resolution(self) -> 'LocusResolution':
        return self.locus.resolution

    def has_resolution(self, resolution: 'LocusResolution') -> bool:
        return self.resolution == resolution

    def _validate(self):
        assert self.component_name is not None and re.match('\w+', self.component_name)


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
    def __init__(self, domain: Optional[str], subdomain: Optional[str], residue: Optional[str]):
        self.domain, self.subdomain, self.residue = domain, subdomain, residue
        self._validate()

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
            raise AssertionError

    def __eq__(self, other: 'Locus') -> bool:
        return self.domain == other.domain and self.subdomain == other.subdomain and self.residue == other.residue

    def _validate(self):
        if self.domain:
            assert re.match("\w+", self.domain)
        if self.subdomain:
            assert re.match("\w+", self.subdomain)
            assert self.domain is not None
        if self.residue:
            assert re.match("\w+", self.residue)

    @property
    def is_empty(self):
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
            raise NotImplementedError


def EmptyLocus():
    return Locus(None, None, None)


class LocusResolution(Enum):
    component = 'component'
    domain    = 'domain'
    subdomain = 'subdomain'
    residue   = 'residue'

    def __lt__(self, other: 'LocusResolution'):
        if self == LocusResolution.component:
            return other != LocusResolution.component
        elif self == LocusResolution.domain:
            return other in [LocusResolution.subdomain, LocusResolution.residue]
        elif self == LocusResolution.subdomain:
            return other == LocusResolution.residue
        else:
            return False


@unique
class SpecSuffix(OrderedEnum):
    mrna    = 'mRNA'
    dna     = 'Gene'
    protein = ''


suffix_to_spec = OrderedDict(
    [
        (SpecSuffix.mrna, MRNASpec),
        (SpecSuffix.dna, GeneSpec),
        (SpecSuffix.protein, ProteinSpec)
    ]
)

spec_to_suffix = OrderedDict((k, v) for v, k in suffix_to_spec.items())


def locus_from_str(locus_str: str) -> Locus:
    def locus_items_from_str(full_locus_str):
        DOMAIN_SUBDOMAIN_RESIDUE_REGEX = '^[\w:-]+\/[\w:-]+\([\w:-]+\)$'
        DOMAIN_RESIDUE_REGEX = '^[\w:-]+\([\w:-]+\)$'
        DOMAIN_SUBDOMAIN_REGEX = '^[\w:-]+\/[\w:-]+$'
        RESIDUE_REGEX = '^\([\w:-]+\)$'
        DOMAIN_REGEX = '^[\w:-]+$'

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
    def spec_from_suffixed_name_and_locus(name: str, struct_index: Optional[int], locus: Locus):
        for suffix in suffix_to_spec:
            if name.endswith(suffix.value):
                name = name[:len(name) - len(suffix.value)]
                return suffix_to_spec[suffix](name, struct_index, locus)
            elif name.lower().endswith(suffix.value.lower()):
                raise SyntaxError('Please use correct capitalization \'{}\' in component \'{}\'.'.format(suffix.value, name))

        raise SyntaxError('Could not parse spec component_name {}'.format(name))

    if not re.match('[A-Za-z][A-Za-z0-9]*(@[\d]+)*(_\[[[A-Za-z0-9]\/\(\)]+\])*', spec_str):
        raise SyntaxError('Spec str {} does not match validating regex.'.format(spec_str))

    DOMAIN_DELIMITER = '_'
    STRUCT_DELIMITER = '@'

    struct_index = None
    items = spec_str.split(DOMAIN_DELIMITER, maxsplit=1)

    if not re.match('[a-zA-Z]', items[0]):
        raise SyntaxError('Spec has to start with letter character.')
    elif STRUCT_DELIMITER in items[0]:
        name, struct_index = items[0].split(STRUCT_DELIMITER)
        struct_index = int(struct_index)
    else:
        name = items[0]

    if len(items) == 1:
        return spec_from_suffixed_name_and_locus(name, struct_index, EmptyLocus())
    elif len(items) == 2:
        return spec_from_suffixed_name_and_locus(name, struct_index, locus_from_str(items[1]))
    else:
        raise SyntaxError('Could not parse spec string {}'.format(spec_str))
