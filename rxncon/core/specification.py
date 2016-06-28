from typecheck import typecheck
import re
from abc import ABCMeta, abstractmethod
from typing import Optional
from enum import Enum


class Specification(metaclass=ABCMeta):
    @typecheck
    def __init__(self, name: str, structure_index: Optional[int], domain: 'Domain'):
        self.name = name
        self.domain = domain
        self.structure_index = structure_index
        self._validate()

    def _validate(self):
        assert self.name is not None and re.match("\w+", self.name)

    def __hash__(self):
        return hash(str(self))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        if self.structure_index:
            return '{0}: {1}@{2}_[{3}/{4}({5})]'.format(str(type(self)), self.structure_index,
                                                        self.name, self.domain, self.subdomain, self.resolution)
        else:
            return '{0}: {1}_[{2}/{3}({4})]'.format(str(type(self)),
                                                    self.name, self.domain, self.subdomain, self.resolution)

    def __eq__(self, other: 'Specification') -> bool:
        return isinstance(other, type(self)) and self.name == other.name \
            and self.domain == other.domain and self.subdomain == other.subdomain and self.residue == other.residue

    def __lt__(self, other: 'Specification') -> bool:
        # None is smaller than something
        if self.name < other.name:
            return True
        elif self.name == other.name:
            return self.domain < other.domain
        return False

    def _validate(self):
        assert self.name is not None and re.match("\w+", self.name)
        if self.domain:
            assert re.match("\w+", self.domain)
        if self.subdomain:
            assert re.match("\w+", self.subdomain)
            assert self.domain is not None
        if self.residue:
            assert re.match("\w+", self.residue)

    @abstractmethod
    def is_equivalent_to(self, other: 'Specification') -> bool:
        pass

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

    def to_domain_resolution(self) -> 'Domain':
        return self.domain

    @property
    def resolution(self) -> 'SpecificationResolution':
        return self.domain.resolution

    def has_resolution(self, resolution: 'SpecificationResolution') -> bool:
        return self.resolution == resolution



class Domain:
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
        return isinstance(other, Domain) and self.domain == other.domain \
               and self.subdomain == other.subdomain and self.residue == other.residue

    @tc.typecheck
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

    @tc.typecheck
    def is_equivalent_to(self, other) -> bool:
        if isinstance(other, Domain) and self.residue and (self.residue == other.residue):
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
        return EmptySpecification(self.name, self.structure_index, Domain(None, None, None))


class ProteinSpecification(Specification):

    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return sfr.string_from_protein_specification(self)

    #@tc.typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, ProteinSpecification) and self.name == other.name \
               and self.domain == other.domain and self.structure_index == other.structure_index

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

    @typecheck
    def is_equivalent_to(self, other: Specification):
        if isinstance(other, ProteinSpecification) and (self.name == other.name) \
                and self.domain.is_equivalent_to(other.domain):
            return True
        else:
            return self == other

    def to_component_specification(self) -> 'ProteinSpecification':
        return ProteinSpecification(self.name, self.structure_index, Domain(None, None, None))


class RnaSpecification(Specification):
    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return sfr.string_from_rna_specification(self)

    #@tc.typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, RnaSpecification) and self.name == other.name \
               and self.domain == other.domain and self.structure_index == other.structure_index

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
                and self.domain.is_equivalent_to(other.domain):
            return True
        else:
            return self == other

    def to_component_specification(self) -> 'RnaSpecification':
        return RnaSpecification(self.name, self.structure_index, Domain(None, None, None))


class DnaSpecification(Specification):
    def __hash__(self):
        return hash(str(self))

    def __str__(self) -> str:
        return sfr.string_from_gene_specification(self)

    #@tc.typecheck
    def __eq__(self, other: Specification) -> bool:
        return isinstance(other, DnaSpecification) and self.name == other.name \
               and self.domain == other.domain and self.structure_index == other.structure_index

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
                and self.domain.is_equivalent_to(other.domain):
            return True
        else:
            return self == other

    def to_component_specification(self) -> 'DnaSpecification':
        return DnaSpecification(self.name, self.structure_index, Domain(None, None, None))


class SpecificationResolution(Enum):
    component = 'component'
    domain    = 'domain'
    subdomain = 'subdomain'
    residue   = 'residue'
