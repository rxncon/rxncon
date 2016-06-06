from typing import Dict
from rxncon.core.specification import Specification
from rxncon.semantics.molecule import MoleculeDefinition, PropertyInstance, AssociationPropertyInstance, \
    create_partner_ass_prop_instance


class Fact:
    pass


class OneParticleFact(Fact):
    def __init__(self, mol_defs: Dict[Specification, MoleculeDefinition], component: Specification,
                 prop_instance: PropertyInstance):
        self.mol_defs = mol_defs
        self.component, self.prop_instance = component, prop_instance

    def __str__(self):
        return 'OneParticleFact:{0}#{1}'.format(str(self.component), str(self.prop_instance))

    def __repr__(self):
        return str(self)

    def complements(self):
        if isinstance(self.prop_instance, AssociationPropertyInstance):
            return self._ass_complements()
        else:
            return [OneParticleFact(self.mol_defs, self.component, x)
                    for x in self.prop_instance.complementary_instances]

    def _ass_complements(self):
        half_bindings = self.prop_instance.complementary_instances
        res = []

        for half_binding in half_bindings:
            partner = create_partner_ass_prop_instance(self.mol_defs, half_binding)
            res.append(TwoParticleFact(self.mol_defs, self.component, half_binding,
                                       half_binding.partner.to_component_specification(), partner))

        return res


class TwoParticleFact(Fact):
    def __init__(self, mol_defs: Dict[Specification, MoleculeDefinition], first_component: Specification,
                 first_prop_instance: PropertyInstance, second_component: Specification,
                 second_prop_instance: PropertyInstance):
        self.mol_defs = mol_defs
        self.first_component, self.first_prop_instance = first_component, first_prop_instance
        self.second_component, self.second_prop_instance = second_component, second_prop_instance
        assert isinstance(self.first_prop_instance, AssociationPropertyInstance)
        assert isinstance(self.second_prop_instance, AssociationPropertyInstance)

    def __str__(self):
        return 'TwoParticleFact:{0}#{1},{2}:{3}'.format(str(self.first_component), str(self.first_prop_instance),
                                                        str(self.second_component), str(self.second_prop_instance))

    def __repr__(self):
        return str(self)

    def complements(self):
        res = []

        for component, instance in [(self.first_component, self.first_prop_instance),
                                    (self.second_component, self.second_prop_instance)]:
            for half_binding in instance.complementary_instances:
                if half_binding.partner:
                    partner = create_partner_ass_prop_instance(self.mol_defs, half_binding)
                    res.append(TwoParticleFact(self.mol_defs, component, half_binding,
                                               half_binding.partner.to_component_specification(), partner))
                else:
                    res.append(OneParticleFact(self.mol_defs, component, half_binding))

        return res



