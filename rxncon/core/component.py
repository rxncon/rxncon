import typecheck as tc
from typing import Optional

import rxncon.syntax.string_from_rxncon as sfr


class Component:
    """Component specifies at some resolution a component of a Reaction or a State."""
    @tc.typecheck
    def __init__(self, name: str, domain: Optional[str], subdomain: Optional[str], residue: Optional[str]):
        self.name = name
        self.domain = domain
        self.subdomain = subdomain
        self.residue = residue

    def __str__(self) -> str:
        return sfr.string_from_component(self)

    @tc.typecheck
    def __eq__(self, other: 'Component') -> bool:
        return self.name == other.name and self.domain == other.domain and self.subdomain == other.subdomain and self.residue == other.residue

    @tc.typecheck
    def is_equivalent_to(self, other: 'Component') -> bool:
        if (self.name == other.name) and self.residue and (self.residue == other.residue):
            return True
        else:
            return self == other

    @tc.typecheck
    def is_subspecification_of(self, other: 'Component') -> bool:
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
    def is_superspecification_of(self, other: 'Component') -> bool:
        if self.is_equivalent_to(other):
            return True

        return other.is_subspecification_of(self)
