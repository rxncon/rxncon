import typecheck as tc
from typing import Optional

import rxncon.syntax.string_from_rxncon as sfr


class Specification:
    @tc.typecheck
    def __init__(self, name: str, domain: Optional[str], subdomain: Optional[str], residue: Optional[str]):
        self.name = name
        self.domain = domain
        self.subdomain = subdomain
        self.residue = residue
        self._validate()

    def _validate(self):
        if self.subdomain:
            assert self.domain is not None

    def __hash__(self) -> int:
        return hash(str(self))

    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return sfr.string_from_component(self)

    @tc.typecheck
    def __eq__(self, other: 'Specification') -> bool:
        return self.name == other.name and self.domain == other.domain and self.subdomain == other.subdomain and self.residue == other.residue

    @tc.typecheck
    def __lt__(self, other: 'Specification'):
        # None is smaller than something
        if self.name < other.name:
            return True
        elif self.name == other.name:
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
    def is_equivalent_to(self, other: 'Specification') -> bool:
        if (self.name == other.name) and self.residue and (self.residue == other.residue):
            return True
        else:
            return self == other

    @tc.typecheck
    def is_subspecification_of(self, other: 'Specification') -> bool:
        if self.is_equivalent_to(other):
            return True

        spec_pairs = zip([self.name, self.domain, self.subdomain, self.residue],
                         [other.name, other.domain, other.subdomain, other.residue])

        for my_property, other_property in spec_pairs:
            if my_property and other_property and my_property != other_property:
                return False

            elif not my_property and other_property:
                return False

        return True

    @tc.typecheck
    def is_superspecification_of(self, other: 'Specification') -> bool:
        if self.is_equivalent_to(other):
            return True

        return other.is_subspecification_of(self)
