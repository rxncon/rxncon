from typecheck import typecheck
from abc import ABCMeta, abstractmethod
from typing import Optional
from enum import Enum
import re

from rxncon.syntax.string_from_rxncon import string_from_domain_information, string_from_protein_specification, \
    string_from_gene_specification, string_from_rna_specification

EMPTY_SPEC = '0'

class Specification(metaclass=ABCMeta):
    @typecheck
    def __init__(self, name: str, structure_index: Optional[int], domain: 'Domain'):
        self.name, self.structure_index, self.domain = name, structure_index, domain
        self._validate()

    def __hash__(self):
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if self.structure_index:
            return '{0}: {1}@{2}_[{3}]'.format(str(type(self)), self.name, self.structure_index, str(self.domain))
        else:
            return '{0}: {1}_[{2}]'.format(str(type(self)), self.name, str(self.domain))

    def __eq__(self, other: 'Specification') -> bool:
        return isinstance(other, type(self)) and self.name == other.name and self.domain == other.domain

    def __lt__(self, other: 'Specification') -> bool:
        if self.name < other.name:
            return True
        elif self.name == other.name:
            return self.domain < other.domain
        return False

    def _validate(self):
        assert self.name is not None and re.match("\w+", self.name)

    def is_subspecification_of(self, other: 'Specification') -> bool:
        if self.is_equivalent_to(other):
            return True

        spec_pairs = zip([self.name, self.domain.domain, self.domain.subdomain, self.domain.residue],
                         [other.name, other.domain.domain, other.domain.subdomain, other.domain.residue])

        for my_property, other_property in spec_pairs:
            if my_property and other_property and my_property != other_property:
                return False

            elif not my_property and other_property:
                return False

        return True

    def is_superspecification_of(self, other: 'Specification') -> bool:
        if self.is_equivalent_to(other):
            return True

        return other.is_subspecification_of(self)

    @abstractmethod
    def to_component_specification(self):
        pass

    def to_dna_component_specification(self) -> 'DnaSpecification':
        return DnaSpecification(self.name, self.structure_index, Domain(None, None, None))

    def to_rna_component_specification(self) -> 'RnaSpecification':
        return RnaSpecification(self.name, self.structure_index, Domain(None, None, None))

    def to_protein_component_specification(self) -> 'ProteinSpecification':
        return ProteinSpecification(self.name, self.structure_index, Domain(None, None, None))

    @property
    def resolution(self) -> 'SpecificationResolution':
        return self.domain.resolution

    def has_resolution(self, resolution: 'SpecificationResolution') -> bool:
        return self.resolution == resolution


class Domain:
    @typecheck
    def __init__(self, domain: Optional[str], subdomain: Optional[str], residue: Optional[str]):
        self.domain, self.subdomain, self.residue = domain, subdomain, residue
        self._validate()

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return string_from_domain_information(self)

    @typecheck
    def __eq__(self, other) -> bool:
        return isinstance(other, Domain) and self.domain == other.domain \
            and self.subdomain == other.subdomain and self.residue == other.residue

    @typecheck
    def __lt__(self, other: 'Domain') -> bool:
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
    def resolution(self) -> 'SpecificationResolution':
        if not self.domain and not self.subdomain and not self.residue:
            return SpecificationResolution.component
        elif self.domain and not self.subdomain and not self.residue:
            return SpecificationResolution.domain
        elif self.domain and self.subdomain and not self.residue:
            return SpecificationResolution.subdomain
        elif self.residue is not None:
            return SpecificationResolution.residue
        else:
            raise NotImplementedError


class EmptySpecification(Specification):
    def _validate(self):
        assert self.name == EMPTY_SPEC

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return EMPTY_SPEC

    def __eq__(self, other: 'Specification') -> bool:
        return isinstance(other, EmptySpecification)

    @typecheck
    def __lt__(self, other: Specification) -> bool:
        if isinstance(other, EmptySpecification):
            return super().__lt__(other)
        elif isinstance(other, ProteinSpecification):
            return True
        elif isinstance(other, RnaSpecification):
            return True
        elif isinstance(other, DnaSpecification):
            return True
        else:
            raise NotImplementedError

    def to_component_specification(self):
        raise NotImplementedError


class ProteinSpecification(Specification):
    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return string_from_protein_specification(self)

    @typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, ProteinSpecification) and self.name == other.name \
            and self.domain == other.domain and self.structure_index == other.structure_index

    @typecheck
    def __lt__(self, other: Specification) -> bool:
        if isinstance(other, ProteinSpecification):
            return super().__lt__(other)
        elif isinstance(other, RnaSpecification):
            return False
        elif isinstance(other, DnaSpecification):
            return False
        elif isinstance(other, EmptySpecification):
            return False
        else:
            raise NotImplementedError

    def to_component_specification(self) -> 'ProteinSpecification':
        return ProteinSpecification(self.name, self.structure_index, Domain(None, None, None))


class RnaSpecification(Specification):
    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return string_from_rna_specification(self)

    @typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, RnaSpecification) and self.name == other.name \
            and self.domain == other.domain and self.structure_index == other.structure_index

    @typecheck
    def __lt__(self, other: Specification) -> bool:
        if isinstance(other, RnaSpecification):
            return super().__lt__(other)
        elif isinstance(other, DnaSpecification):
            return False
        elif isinstance(other, ProteinSpecification):
            return True
        elif isinstance(other, EmptySpecification):
            return False
        else:
            raise NotImplementedError

    def to_component_specification(self) -> 'RnaSpecification':
        return RnaSpecification(self.name, self.structure_index, Domain(None, None, None))


class DnaSpecification(Specification):
    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return string_from_gene_specification(self)

    @typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, DnaSpecification) and self.name == other.name \
               and self.domain == other.domain and self.structure_index == other.structure_index

    @typecheck
    def __lt__(self, other: Specification):
        if isinstance(other, DnaSpecification):
            return super().__lt__(other)
        elif isinstance(other, RnaSpecification):
            return True
        elif isinstance(other, ProteinSpecification):
            return True
        elif isinstance(other, EmptySpecification):
            return False
        else:
            raise NotImplementedError

    def to_component_specification(self) -> 'DnaSpecification':
        return DnaSpecification(self.name, self.structure_index, Domain(None, None, None))


class SpecificationResolution(Enum):
    component = 'component'
    domain    = 'domain'
    subdomain = 'subdomain'
    residue   = 'residue'
