from typing import Dict, List, Optional

from rxncon.core.specification import Specification
from rxncon.core.contingency import ContingencyType
from rxncon.core.effector import Effector, StateEffector, AndEffector, OrEffector, NotEffector
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.core.reaction import Reactant as RxnConReactant, Reaction
from rxncon.semantics.molecule import MoleculeDefinition, MoleculeInstance, MutualExclusivityError
from rxncon.semantics.molecule_from_rxncon import mol_defs_from_rxncon_sys
from rxncon.semantics.elemental import elemental_from_state, OneParticleElemental, TwoParticleElemental, Elemental
from rxncon.venntastic.sets import PropertySet, Intersection, Union, Complement, UniversalSet, gram_schmidt_disjunctify
from rxncon.simulation.rule_based.rule_based_model import Rule, Complex, Arrow, Parameter, Binding
from rxncon.util.utils import transform_set_expression


def state_set_from_effector(effector: Effector):
    if isinstance(effector, StateEffector):
        return PropertySet(effector.expr)
    elif isinstance(effector, NotEffector):
        return Complement(state_set_from_effector(effector.expr))
    elif isinstance(effector, AndEffector):
        return Intersection(state_set_from_effector(effector.left_expr),
                            state_set_from_effector(effector.right_expr))
    elif isinstance(effector, OrEffector):
        return Union(state_set_from_effector(effector.left_expr),
                     state_set_from_effector(effector.right_expr))
    else:
        raise NotImplementedError


def elemental_set_from_state_set(mol_defs: Dict[Specification, MoleculeDefinition], state_set):
    return transform_set_expression(state_set, lambda x: elemental_from_state(mol_defs, x))


def rules_from_reaction(rxnconsys: RxnConSystem, reaction: Reaction) -> List[Rule]:
    mol_defs = mol_defs_from_rxncon_sys(rxnconsys)
    state_set = UniversalSet()
    for contingency in rxnconsys.strict_contingencies_for_reaction(reaction):
        if contingency.type == ContingencyType.requirement:
            state_set = Intersection(state_set, state_set_from_effector(contingency.effector))
        elif contingency.type == ContingencyType.inhibition:
            state_set = Intersection(state_set, Complement(state_set_from_effector(contingency.effector)))
        else:
            raise NotImplementedError

    # Convert to set of elementals
    elemental_set = elemental_set_from_state_set(mol_defs, state_set)

    # Simplify the complements
    elemental_set = elemental_set.simplified_form()

    # Get the different solutions to the contingencies, still overlapping
    overlapping_solns = elemental_set.to_union_list_form()

    # Disjunctify them
    disjunct_solns = gram_schmidt_disjunctify(overlapping_solns)

    # Simplify once more to root out the last complements
    disjunct_solns = [disjunct_soln.simplified_form() for disjunct_soln in disjunct_solns]

    rules = []

    for soln in disjunct_solns:
        possible_rule = rule_from_reaction_and_contingency_soln(mol_defs, reaction, soln)
        if possible_rule:
            rules.append(possible_rule)

    return rules


def rule_from_reaction_and_contingency_soln(mol_defs, reaction, soln) -> Optional[Rule]:
    assert len(soln.to_nested_list_form()) == 1

    lhs_reacting_elementals = [elemental_from_state(mol_defs, state)
                               for reactant in reaction.reactants_pre for state in reactant.states]

    rhs_reacting_elementals = [elemental_from_state(mol_defs, state)
                               for reactant in reaction.reactants_post for state in reactant.states]

    lhs_reactants = complexes_from_elementals(mol_defs, lhs_reacting_elementals, soln.to_nested_list_form()[0])
    rhs_reactants = complexes_from_elementals(mol_defs, rhs_reacting_elementals, soln.to_nested_list_form()[0])

    arrow_type = Arrow.irreversible
    parameter = Parameter('x', None)

    if lhs_reactants and rhs_reactants:
        return Rule(set(lhs_reactants), set(rhs_reactants), arrow_type, set(parameter))
    else:
        return None


def complexes_from_elementals(mol_defs, reacting_elementals: List[Elemental], background_elementals: List[Elemental]):
    class MoleculeDict(dict):
        def __missing__(self, key):
            self[key] = MoleculeInstance(mol_defs[key], set(), set(), None)
            return self[key]

    molecules = MoleculeDict()
    bindings = []
    reacting_components = []

    for reacting_elemental in reacting_elementals:
        if isinstance(reacting_elemental, OneParticleElemental):
            reacting_components.append(reacting_elemental.component)
            molecules[reacting_elemental.component].add_property(reacting_elemental.prop_instance)
        elif isinstance(reacting_elemental, TwoParticleElemental):
            reacting_components += [reacting_elemental.first_component, reacting_elemental.second_component]
            molecules[reacting_elemental.first_component].add_property(reacting_elemental.first_prop_instance)
            molecules[reacting_elemental.second_component].add_property(reacting_elemental.second_prop_instance)

    try:
        for elem in [x for x in background_elementals if isinstance(x, TwoParticleElemental)]:  # type: TwoParticleElemental
            molecules[elem.first_component].add_property(elem.first_prop_instance)
            molecules[elem.second_component].add_property(elem.second_prop_instance)
            bindings.append(Binding(elem.first_prop_instance.property_def.spec,
                                    elem.second_prop_instance.property_def.spec))

        for elem in [x for x in background_elementals if isinstance(x, OneParticleElemental)]:  # type: OneParticleElemental
            molecules[elem.component].add_property(elem.prop_instance)

    except MutualExclusivityError:
        return None

    complexes = []

    for reacting_component in reacting_components:
        if reacting_component not in molecules.keys():
            # This can only be the case if the reacting component has already been joined in a complex with
            # another reacting component, such as in the reaction A + B --> A.B
            continue

        the_components = [reacting_component]
        the_molecules = [molecules[reacting_component]]
        the_bindings = []

        continue_trying = True
        while continue_trying:
            continue_trying = False

            for binding in bindings:
                if any(binding.left_partner.is_subspecification_of(x) for x in the_components):
                    continue_trying = True
                    the_bindings.append(binding)
                    bindings.remove(binding)
                    if not binding.right_partner.to_component_specification() in the_components:
                        the_components.append(binding.right_partner.to_component_specification())
                        the_molecules.append(molecules[binding.right_partner.to_component_specification()])
                        molecules.pop(binding.right_partner.to_component_specification())

                elif any(binding.right_partner.is_subspecification_of(x) for x in the_components):
                    continue_trying = True
                    the_bindings.append(binding)
                    bindings.remove(binding)
                    if not binding.left_partner.to_component_specification() in the_components:
                        the_components.append(binding.left_partner.to_component_specification())
                        the_molecules.append(molecules[binding.left_partner.to_component_specification()])
                        molecules.pop(binding.left_partner.to_component_specification())

        complexes.append(Complex(set(the_molecules), set(the_bindings)))

    return complexes








