import typecheck as tc
import re
from abc import ABCMeta, abstractmethod
from typing import Optional
from enum import Enum

import rxncon.syntax.string_from_rxncon as sfr


class Specification(metaclass=ABCMeta):

    @tc.typecheck
    def __init__(self, name: str, structure_index: Optional[int], spec_resolution: 'DomainDefinition'):
        self.name = name
        self.spec_resolution = spec_resolution
        self.structure_index = structure_index
        self._validate()

    def _validate(self):
        assert self.name is not None and re.match("\w+", self.name)

    @abstractmethod
    def __hash__(self) -> int:
        pass

    def __repr__(self) -> str:
        return str(self)

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def __eq__(self, other: 'Specification') -> bool:
        pass

    def __lt__(self, other: 'Specification') -> bool:
        # None is smaller than something
        if self.name < other.name:
            return True
        elif self.name == other.name:
            return self.spec_resolution < other.spec_resolution
        return False

    @abstractmethod
    def is_equivalent_to(self, other: 'Specification') -> bool:
        pass

    def is_subspecification_of(self, other: 'Specification') -> bool:
        if self.is_equivalent_to(other):
            return True

        spec_pairs = zip([self.name, self.spec_resolution.domain, self.spec_resolution.subdomain, self.spec_resolution.residue],
                         [other.name, other.spec_resolution.domain, other.spec_resolution.subdomain, other.spec_resolution.residue])

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
        return DnaSpecification(self.name, self.structure_index, DomainDefinition(None, None, None))

    def to_rna_component_specification(self) -> 'RnaSpecification':
        return RnaSpecification(self.name, self.structure_index, DomainDefinition(None, None, None))

    def to_protein_component_specification(self) -> 'ProteinSpecification':
        return ProteinSpecification(self.name, self.structure_index, DomainDefinition(None, None, None))

    def to_domain_resolution(self) -> 'DomainDefinition':
        return self.spec_resolution

    @property
    def resolution(self) -> 'SpecificationResolution':
        return self.spec_resolution.resolution

    def has_resolution(self, resolution: 'SpecificationResolution') -> bool:
        return self.resolution == resolution



class DomainDefinition:
    @tc.typecheck
    def __init__(self, domain: Optional[str], subdomain: Optional[str], residue: Optional[str]):
        self.domain = domain
        self.subdomain = subdomain
        self.residue = residue
        self._validate()

    def _validate(self):
        if self.domain:
            assert re.match("\w+", self.domain)
        if self.subdomain:
            assert re.match("\w+", self.subdomain)
            assert self.domain is not None
        if self.residue:
            assert re.match("\w+", self.residue)

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return sfr.string_from_domain_information(self)

    @tc.typecheck
    def __eq__(self, other) -> bool:
        return isinstance(other, DomainDefinition) and self.domain == other.domain \
               and self.subdomain == other.subdomain and self.residue == other.residue

    @tc.typecheck
    def __lt__(self, other: 'DomainDefinition') -> bool:
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

    @tc.typecheck
    def is_equivalent_to(self, other) -> bool:
        if isinstance(other, DomainDefinition) and self.residue and (self.residue == other.residue):
            return True
        else:
            return self == other

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
        assert self.name is not None and re.match("0", self.name)

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        return sfr.string_from_protein_specification(self)

    def __eq__(self, other: 'Specification') -> bool:
        return isinstance(other, EmptySpecification) and self.name == other.name

    @tc.typecheck
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

    @tc.typecheck
    def is_equivalent_to(self, other: Specification) -> bool:
        return self == other

    def to_component_specification(self) -> 'EmptySpecification':
        return EmptySpecification(self.name, self.structure_index, DomainDefinition(None, None, None))


class ProteinSpecification(Specification):

    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return sfr.string_from_protein_specification(self)

    #@tc.typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, ProteinSpecification) and self.name == other.name \
               and self.spec_resolution == other.spec_resolution and self.structure_index == other.structure_index

    @tc.typecheck
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

    @tc.typecheck
    def is_equivalent_to(self, other: Specification) -> bool:
        if isinstance(other, ProteinSpecification) and (self.name == other.name) \
                and self.spec_resolution.is_equivalent_to(other.spec_resolution):
            return True
        else:
            return self == other

    def to_component_specification(self) -> 'ProteinSpecification':
        return ProteinSpecification(self.name, self.structure_index, DomainDefinition(None, None, None))


class RnaSpecification(Specification):

    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return sfr.string_from_rna_specification(self)

    #@tc.typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, RnaSpecification) and self.name == other.name \
               and self.spec_resolution == other.spec_resolution and self.structure_index == other.structure_index

    @tc.typecheck
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

    @tc.typecheck
    def is_equivalent_to(self, other: Specification) -> bool:
        if isinstance(other, RnaSpecification) and (self.name == other.name) \
                and self.spec_resolution.is_equivalent_to(other.spec_resolution):
            return True
        else:
            return self == other

    def to_component_specification(self) -> 'RnaSpecification':
        return RnaSpecification(self.name, self.structure_index, DomainDefinition(None, None, None))


class DnaSpecification(Specification):

    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return sfr.string_from_gene_specification(self)

    #@tc.typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, DnaSpecification) and self.name == other.name \
               and self.spec_resolution == other.spec_resolution and self.structure_index == other.structure_index

    @tc.typecheck
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

    @tc.typecheck
    def is_equivalent_to(self, other: Specification) -> bool:
        if isinstance(other, DnaSpecification) and (self.name == other.name) \
                and self.spec_resolution.is_equivalent_to(other.spec_resolution):
            return True
        else:
            return self == other

    def to_component_specification(self) -> 'DnaSpecification':
        return DnaSpecification(self.name, self.structure_index, DomainDefinition(None, None, None))


class SpecificationResolution(Enum):
    component = 'component'
    domain    = 'domain'
    subdomain = 'subdomain'
    residue   = 'residue'
