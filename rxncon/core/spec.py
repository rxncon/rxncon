from collections import OrderedDict
from typecheck import typecheck
from abc import ABCMeta, abstractmethod
from typing import Optional, Union
from enum import Enum, unique
import re

from rxncon.util.utils import OrderedEnum

EMPTY_MOL_SPEC = '0'

class Spec(metaclass=ABCMeta):
    pass


class MolSpec(Spec, metaclass=ABCMeta):
    @typecheck
    def __init__(self, component_name: str, struct_index: Optional[int], locus: 'Locus'):
        self.component_name, self.struct_index, self.locus = component_name, struct_index, locus
        self._validate()

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        @typecheck
        def _string_from_spec(spec: MolSpec) -> str:
            def struct_name(spec: MolSpec, suffix: 'SpecSuffix'):
                if spec.struct_index:
                    return "{0}{1}@{2}".format(spec.component_name, suffix.value, spec.struct_index)
                else:
                    return "{0}{1}".format(spec.component_name, suffix.value)

            suffix = spec_to_suffix[type(spec)]

            if str(spec.locus):
                return '{0}_[{1}]'.format(struct_name(spec, suffix), str(spec.locus))
            else:
                return '{0}'.format(struct_name(spec, suffix))

        return _string_from_spec(self)

    @typecheck
    def __eq__(self, other: 'MolSpec') -> bool:
        return isinstance(other, type(self)) and self.component_name == other.component_name and self.locus == other.locus \
            and self.struct_index == other.struct_index

    @typecheck
    def __lt__(self, other: 'MolSpec') -> bool:
        if self.component_name < other.component_name:
            return True
        elif self.component_name == other.component_name and self.locus == other.locus:
            return self.struct_index < other.struct_index
        elif self.component_name == other.component_name:
            return self.locus < other.locus
        else:
            return False

    def _validate(self):
        assert self.component_name is not None and re.match('\w+', self.component_name)

    @typecheck
    def is_equivalent_to(self, other: 'MolSpec') -> bool:
        return self == other or type(self) == type(other) and self.locus.residue == other.locus.residue and \
            self.struct_index == other.struct_index

    @typecheck
    def is_subspec_of(self, other: 'MolSpec') -> bool:
        if self.is_equivalent_to(other):
            return True

        spec_pairs = zip([self.component_name, self.locus.domain, self.locus.subdomain, self.locus.residue],
                         [other.component_name, other.locus.domain, other.locus.subdomain, other.locus.residue])

        for my_property, other_property in spec_pairs:
            if my_property and other_property and my_property != other_property:
                return False

            elif not my_property and other_property:
                return False

        return True

    @typecheck
    def is_superspec_of(self, other: 'MolSpec') -> bool:
        if self.is_equivalent_to(other):
            return True

        return other.is_subspec_of(self)

    @property
    @typecheck
    def is_component_spec(self) -> bool:
        return self.has_resolution(LocusResolution.component)

    @abstractmethod
    @typecheck
    def to_component_spec(self) -> 'MolSpec':
        pass

    @typecheck
    def to_protein_component_spec(self) -> 'ProteinSpec':
        return ProteinSpec(self.component_name, self.struct_index, EmptyLocus())

    @typecheck
    def to_dna_component_spec(self) -> 'DnaSpec':
        return DnaSpec(self.component_name, self.struct_index, EmptyLocus())

    @typecheck
    def to_mrna_component_spec(self) -> 'MRnaSpec':
        return MRnaSpec(self.component_name, self.struct_index, EmptyLocus())

    @property
    @typecheck
    def resolution(self) -> 'LocusResolution':
        return self.locus.resolution

    @typecheck
    def has_resolution(self, resolution: 'LocusResolution') -> bool:
        return self.resolution == resolution


class EmptyMolSpec(MolSpec):
    def __init__(self):
        super().__init__(EMPTY_MOL_SPEC, None, EmptyLocus())

    def _validate(self):
        assert self.component_name == EMPTY_MOL_SPEC
        assert self.locus.is_empty

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return EMPTY_MOL_SPEC

    @typecheck
    def __eq__(self, other: 'MolSpec') -> bool:
        return isinstance(other, EmptyMolSpec)

    @typecheck
    def __lt__(self, other: MolSpec) -> bool:
        if isinstance(other, EmptyMolSpec):
            return False
        elif isinstance(other, ProteinSpec):
            return True
        elif isinstance(other, MRnaSpec):
            return True
        elif isinstance(other, DnaSpec):
            return True
        else:
            raise NotImplementedError

    def to_component_spec(self):
        raise AssertionError

    def has_resolution(self, resolution: 'LocusResolution'):
        return True


class ProteinSpec(MolSpec):
    def __hash__(self) -> int:
        return hash(str(self))

    @typecheck
    def __eq__(self, other: MolSpec) -> bool:
        return isinstance(other, ProteinSpec) and self.component_name == other.component_name \
            and self.locus == other.locus and self.struct_index == other.struct_index

    @typecheck
    def __lt__(self, other: MolSpec) -> bool:
        if isinstance(other, ProteinSpec):
            return super().__lt__(other)
        elif isinstance(other, MRnaSpec):
            return False
        elif isinstance(other, DnaSpec):
            return False
        elif isinstance(other, EmptyMolSpec):
            return False
        else:
            raise NotImplementedError

    def to_component_spec(self) -> 'ProteinSpec':
        return ProteinSpec(self.component_name, self.struct_index, EmptyLocus())


class MRnaSpec(MolSpec):
    def __hash__(self) -> int:
        return hash(str(self))

    @typecheck
    def __eq__(self, other: MolSpec) -> bool:
        return isinstance(other, MRnaSpec) and self.component_name == other.component_name \
            and self.locus == other.locus and self.struct_index == other.struct_index

    @typecheck
    def __lt__(self, other: MolSpec) -> bool:
        if isinstance(other, MRnaSpec):
            return super().__lt__(other)
        elif isinstance(other, DnaSpec):
            return False
        elif isinstance(other, ProteinSpec):
            return True
        elif isinstance(other, EmptyMolSpec):
            return False
        else:
            raise NotImplementedError

    def to_component_spec(self) -> 'MRnaSpec':
        return MRnaSpec(self.component_name, self.struct_index, EmptyLocus())


class DnaSpec(MolSpec):
    def __hash__(self) -> int:
        return hash(str(self))

    @typecheck
    def __eq__(self, other: MolSpec) -> bool:
        return isinstance(other, DnaSpec) and self.component_name == other.component_name \
            and self.locus == other.locus and self.struct_index == other.struct_index

    @typecheck
    def __lt__(self, other: MolSpec):
        if isinstance(other, DnaSpec):
            return super().__lt__(other)
        elif isinstance(other, MRnaSpec):
            return True
        elif isinstance(other, ProteinSpec):
            return True
        elif isinstance(other, EmptyMolSpec):
            return False
        else:
            raise NotImplementedError

    def to_component_spec(self) -> 'DnaSpec':
        return DnaSpec(self.component_name, self.struct_index, EmptyLocus())


class Locus:
    @typecheck
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

    @typecheck
    def __eq__(self, other: 'Locus') -> bool:
        return self.domain == other.domain and self.subdomain == other.subdomain and self.residue == other.residue

    @typecheck
    def __lt__(self, other: 'Locus') -> bool:
        if self.domain is None and other.domain is not None:
            return True
        if other.domain is not None and other.domain is not None \
                and self.domain < other.domain:
            return True
        if self.subdomain is None and other.subdomain is not None:
            return True
        if self.subdomain is not None and other.subdomain is not None \
                and self.subdomain < other.subdomain:
            return True
        if self.residue is None and other.residue is not None:
            return True
        if self.residue is not None and other.residue is not None \
                and self.residue < other.residue:
            return True
        return False

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

    @typecheck
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
    dna     = 'DNA'
    protein = ''


suffix_to_spec = OrderedDict(
    [
        (SpecSuffix.mrna, MRnaSpec),
        (SpecSuffix.dna, DnaSpec),
        (SpecSuffix.protein, ProteinSpec)
    ]
)

spec_to_suffix = OrderedDict((k, v) for v, k in suffix_to_spec.items())


@typecheck
def locus_from_string(locus_str: str) -> Locus:
    def locus_items_from_string(full_locus_str):
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

    return Locus(*locus_items_from_string(locus_str.strip('[]')))


class BondSpec(Spec):
    @typecheck
    def __init__(self, first: MolSpec, second: MolSpec):
        self.first, self.second = sorted([first, second])

    @typecheck
    def __eq__(self, other: 'BondSpec') -> bool:
        return self.first == other.first and self.second == other.second

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return 'BondSpec<{0}, {1}>'.format(str(self.first), str(self.second))

    def __repr__(self) -> str:
        return str(self)


@typecheck
def mol_spec_from_string(spec_str: str) -> MolSpec:
    @typecheck
    def spec_from_suffixed_name_and_locus(name: str, struct_index: Optional[int], locus: Locus):
        for suffix in suffix_to_spec:
            if name.endswith(suffix.value):
                name = name[:len(name) - len(suffix.value)]
                return suffix_to_spec[suffix](name, struct_index, locus)

        raise AssertionError('Could not parse spec component_name {}'.format(name))

    DOMAIN_DELIMITER = '_'
    STRUCT_DELIMITER = '@'

    struct_index = None
    items = spec_str.split(DOMAIN_DELIMITER, maxsplit=1)

    if items[0].startswith(EMPTY_MOL_SPEC):
        if items[0] != EMPTY_MOL_SPEC:
            raise SyntaxError('Only the EmptySpec can start with {}'.format(EMPTY_MOL_SPEC))
        return EmptyMolSpec()
    elif STRUCT_DELIMITER in items[0]:
        name, struct_index = items[0].split(STRUCT_DELIMITER)
        struct_index = int(struct_index)
    else:
        name = items[0]

    if len(items) == 1:
        return spec_from_suffixed_name_and_locus(name, struct_index, EmptyLocus())
    elif len(items) == 2:
        return spec_from_suffixed_name_and_locus(name, struct_index, locus_from_string(items[1]))
    else:
        raise SyntaxError('Could not parse spec string {}'.format(spec_str))


@typecheck
def bond_spec_from_string(spec_str: str) -> Spec:
    first, second = sorted([mol_spec_from_string(x) for x in spec_str.split('~')])

    if isinstance(first, EmptyMolSpec):
        return second
    else:
        return BondSpec(first, second)


@typecheck
def spec_from_string(spec_str: str) -> Spec:
    return bond_spec_from_string(spec_str) if '~' in spec_str else mol_spec_from_string(spec_str)

