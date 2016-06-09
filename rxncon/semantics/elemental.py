from typing import Dict, Optional, List
from rxncon.core.specification import Specification
from rxncon.semantics.molecule import MoleculeDefinition, PropertyInstance, AssociationPropertyInstance, \
    create_partner_ass_prop_instance, ModificationPropertyInstance
from rxncon.core.state import State, CovalentModificationState, InteractionState
from rxncon.simulation.rule_based.molecule_from_string import property_ins_from_string
from rxncon.core.reaction import Reactant


EMPTY = '0'

class Elemental:
    pass


class OneParticleElemental(Elemental):
    def __init__(self, mol_defs: Dict[Specification, MoleculeDefinition], component: Specification,
                 prop_instance: Optional[PropertyInstance]):
        self.mol_defs = mol_defs
        self.component, self.prop_instance = component, prop_instance

    def __str__(self):
        return 'OneParticleElemental:{0}#{1}'.format(str(self.component), str(self.prop_instance))

    def __repr__(self):
        return str(self)

    def complements(self):
        if not self.prop_instance:
            return [self]
        elif isinstance(self.prop_instance, AssociationPropertyInstance):
            return self._ass_complements()
        else:
            return [OneParticleElemental(self.mol_defs, self.component, x)
                    for x in self.prop_instance.complementary_instances]

    def _ass_complements(self):
        half_bindings = self.prop_instance.complementary_instances
        res = []

        for half_binding in half_bindings:
            partner = create_partner_ass_prop_instance(self.mol_defs, half_binding)
            res.append(TwoParticleElemental(self.mol_defs, self.component, half_binding,
                                            half_binding.partner.to_component_specification(), partner))

        return res


class TwoParticleElemental(Elemental):
    def __init__(self, mol_defs: Dict[Specification, MoleculeDefinition], first_component: Specification,
                 first_prop_instance: PropertyInstance, second_component: Specification,
                 second_prop_instance: PropertyInstance):
        self.mol_defs = mol_defs
        self.first_component, self.first_prop_instance = first_component, first_prop_instance
        self.second_component, self.second_prop_instance = second_component, second_prop_instance
        assert isinstance(self.first_prop_instance, AssociationPropertyInstance)
        assert isinstance(self.second_prop_instance, AssociationPropertyInstance)

    def __str__(self):
        return 'TwoParticleElemental:{0}#{1},{2}:{3}'.format(str(self.first_component), str(self.first_prop_instance),
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
                    res.append(TwoParticleElemental(self.mol_defs, component, half_binding,
                                                    half_binding.partner.to_component_specification(), partner))
                else:
                    res.append(OneParticleElemental(self.mol_defs, component, half_binding))

        return res


def elemental_from_state(mol_defs: Dict[Specification, MoleculeDefinition], state: State) -> Elemental:
    if isinstance(state, CovalentModificationState):
        component = state.substrate.to_component_specification()
        return OneParticleElemental(mol_defs, component,
                                    property_ins_from_string(mol_defs[component],
                                                             'mod/{0}:{1}'.format(str(state.substrate),
                                                                                  state.modifier.value)))
    elif isinstance(state, InteractionState) and \
            not state.first_component.name == EMPTY and not state.second_component.name == EMPTY:
        first_component = state.first_component.to_component_specification()
        second_component = state.second_component.to_component_specification()
        return TwoParticleElemental(mol_defs, first_component,
                                    property_ins_from_string(mol_defs[first_component],
                                                             'ass/{0}:{1}'.format(str(state.first_component),
                                                                                  str(state.second_component))),
                                    second_component,
                                    property_ins_from_string(mol_defs[second_component],
                                                             'ass/{0}:{1}'.format(str(state.second_component),
                                                                                  str(state.first_component))))
    elif isinstance(state, InteractionState):
        if state.first_component.name == EMPTY:
            the_component = state.second_component
        else:
            the_component = state.first_component

        return OneParticleElemental(mol_defs, the_component.to_component_specification(),
                                    property_ins_from_string(mol_defs[the_component.to_component_specification()],
                                                             'ass/{0}:'.format(str(the_component))))


    else:
        raise AssertionError

def elementals_from_reactant(mol_defs: Dict[Specification, MoleculeDefinition], reactant: Reactant) -> List[Elemental]:
    if reactant.states:
        return [elemental_from_state(mol_defs, x) for x in reactant.states]
    else:
        return [OneParticleElemental(mol_defs, reactant.component, None)]
