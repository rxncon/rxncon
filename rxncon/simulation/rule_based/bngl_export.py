from typing import List, Optional

class Unmodified:
    unmodified = "U"

class DefaultDomain:
    """
    we have to define default domains in BNGL
    for Association reactions this will be Assoc[Partner]
    for relocalisation this will be loc
    for modifications this will be Mod[Modifier]
    """
    default_interacton_domain = "Assoc"
    default_modification_domain = "Mod"
    default_localisation_domain = "loc"


class BNGLSystem:
    def __init__(self, rule_based_system: List[Rule]):
        self.rule_based_system = rule_based_system

        self.molecules = []  # type: List[Molecule]
        self.parameters = []  # type: List[Parameter]


    def bngl_from_rule_based_model(self):
        rule_section = RuleSection(self.rule_based_system)
        bngl = rule_section._create_rule_section()
        self.molecules = bngl.molecules
        self.parameters = bngl.parameters

    def write_to_file(self):
        self._write_molecule_type_section()
        self._write_seeded_species_section()
        self._write_parameter_section()
        self.rule_section.write_to_file()

    def _write_molecule_type_section(self):
        """

        In this section all domains and modifications are considered
        Ste7(ALS359~U~P, AssocSte11)

        The names of the different molecules of the reactant should be sorted
        A().B()
        """
        pass

    def _write_seeded_species_section(self):
        """
        all species should be initilized as not modified
        Ste7(ALS359~U)
        """
        pass

    def _write_parameter_section(self):
        """
        all parameters of the system
        """
        pass

    def has_molecule(self, other: str):
        for mol in self.moleculs:
            if mol.name == other:
                return True
        return False


class Molecule:
    def __init__(self):
        self.name = None
        self.occupied_binding_domain = []
        self.free_binding_domain = []
        self.occupied_modification_domain = []
        self.free_modification_domain = []
        self.relocalisation = []


class Parameter:
    def __init__(self, forward_name: Optional[str], reverse_name: Optional[str]):
        self.forward_name = forward_name
        self.reverse_name = reverse_name
        self.forward_value = 1
        self.reverse_value = 1


class RuleSection:

    def __init__(self, rules: List[Rule], molecules: List[Molecule], parmeters: List[Parameter]):
        self.rules = rules  # type: List[Rule]
        self.molecules = []  # type: List[Molecule]
        self.parameters = []  # type: List[Parameter]

        self._create_rule_section()

    def _create_rule_section(self):
        pass

    def source_to_str(self):
        pass

    def product_to_str(self):
        pass

    def reactant_to_str(self):
        """
        # 1. Localisation domains
        # 2. modification domains
        # 3. Binding domains

        A(1,2,3) + B(1,2,3) <-> A(1,2,3!0).B(1,2,3!0)

        """
        pass

    @property
    def arrow(self):
        """Returns arrow depending on rule reversibility."""
        if self.rule.directionality == 2:  # irreversible
            return '->'
        elif self.rule.directionality == 1:  # reversible
            return '<->'
