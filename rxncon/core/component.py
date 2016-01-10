from typing import Optional

import rxncon.syntax.string_from_rxncon as sfr

class Component:
    def __init__(self, name: str, domain: Optional[str], subdomain: Optional[str], residue: Optional[str]):
        self.name = name
        self.domain = domain
        self.subdomain = subdomain
        self.residue = residue

    def __str__(self) -> str:
        return sfr.string_from_component(self)

    def __eq__(self, other: 'Component') -> bool:
        assert isinstance(other, Component)
        return self.name == other.name and self.domain == other.domain and self.subdomain == other.subdomain and self.residue == other.residue

    def is_equivalent_to(self, other: 'Component'):
        assert isinstance(other, Component)

        if self.residue and self.residue == other.residue:
            return True

        else:
            return self == other

    def is_subspecification_of(self, other: 'Component') -> bool:
        assert isinstance(other, Component)

        if self.is_equivalent_to(other):
            return True

        # @todo Implement this.

    def is_superspecification_of(self, other: 'Component') -> bool:
        assert isinstance(other, Component)

        if self.is_equivalent_to(other):
            return True

        # @todo Implement this.
