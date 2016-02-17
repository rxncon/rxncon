import copy
import typing as tg
from collections import defaultdict

import rxncon.semantics.molecule as mol
import rxncon.core.rxncon_system as rxs
import rxncon.core.state as sta
import rxncon.core.reaction as rxn
import rxncon.core.specification as spe
import rxncon.core.contingency as con
import rxncon.venntastic.sets as venn
import rxncon.core.effector as eff


# CREATION OF MOLECULE DEFINITIONS FROM REACTIONS
class MoleculeDefinitionSupervisor:
    def __init__(self, rxnconsys: rxs.RxnConSystem):
        self.rxnconsys = rxnconsys
        self.molecule_definitions = {}
        self._generate_molecule_definitions()
        self.molecules = self.molecule_definitions.keys()

    def molecule_definition_for_name(self, name: str) -> mol.MoleculeDefinition:
        return self.molecule_definitions[name]

    def _generate_molecule_definitions(self):
        names = set()

        name_to_assoc_defs = defaultdict(set)
        name_to_mod_defs = defaultdict(set)
        name_to_locs = defaultdict(set)

        for reaction in self.rxnconsys.reactions:
            for state in [x for x in [reaction.source, reaction.product] if x]:
                if isinstance(state, sta.CovalentModificationState):
                    names.add(state.substrate.name)
                    names.add(reaction.subject.name)

                    mod_def = _mod_def_from_state_and_reaction(state, reaction)
                    _update_defs(name_to_mod_defs[state.substrate.name], mod_def)

                elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
                    names.add(state.first_component.name)
                    names.add(state.second_component.name)

                    assoc_defs = _assoc_defs_from_state(state)
                    _update_defs(name_to_assoc_defs[state.first_component.name], assoc_defs[0])
                    _update_defs(name_to_assoc_defs[state.second_component.name], assoc_defs[1])

                elif isinstance(state, sta.TranslocationState):
                    names.add(state.substrate.name)
                    names.add(reaction.subject.name)
                    name_to_locs[state.substrate.name].add(state.compartment)

                else:
                    raise NotImplementedError

        for name in names:
            universal_specification = spe.Specification(name, None, None, None)
            self.molecule_definitions[name] = mol.MoleculeDefinition(universal_specification,
                                                                     name_to_mod_defs[name],
                                                                     name_to_assoc_defs[name],
                                                                     mol.LocalizationDefinition(name_to_locs[name]))


# MATCHING MOLECULE DEFINITIONS WITH SETS OF STATES INTO ASSOC/MOD/LOC INSTANCES (_NOT_ MOLECULE INSTANCES YET)
def set_of_instances_from_molecule_def_and_set_of_states(mol_def: mol.MoleculeDefinition, set_of_states: venn.Set) -> venn.Set:
    if set_of_states.is_equivalent_to(venn.EmptySet()):
        raise NotImplementedError

    elif set_of_states.is_equivalent_to(venn.UniversalSet()):
        return venn.UniversalSet()

    elif isinstance(set_of_states, venn.PropertySet):
        instances = _instances(mol_def, set_of_states.value, negate=False)
        if not instances:
            return venn.UniversalSet()
        else:
            return venn.nested_expression_from_list_and_binary_op([venn.PropertySet(x) for x in instances], venn.Union)

    elif isinstance(set_of_states, venn.Complement):
        assert isinstance(set_of_states.expr, venn.PropertySet)
        instances = _instances(mol_def, set_of_states.expr.value, negate=True)
        if not instances:
            return venn.UniversalSet()
        else:
            return venn.nested_expression_from_list_and_binary_op([venn.PropertySet(x) for x in instances], venn.Union)

    elif isinstance(set_of_states, venn.Union):
        return venn.Union(
            set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.left_expr),
            set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.right_expr)
        )

    elif isinstance(set_of_states, venn.Intersection):
        return venn.Intersection(
            set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.left_expr),
            set_of_instances_from_molecule_def_and_set_of_states(mol_def, set_of_states.right_expr)
        )


# CREATING A CONCRETE MOLECULE INSTANCE FROM THE SET OF ASSOC/MOD/LOC INSTANCES
def molecule_instance_from_molecule_def_and_set_of_instances(mol_def: mol.MoleculeDefinition, set_of_instances: venn.Set) -> mol.ModificationInstance:
    instances = set_of_instances.to_nested_list_form()
    assert len(instances) == 1

    instances = instances[0]
    assoc_instances = set()
    mod_instances = set()
    loc = None

    for instance in instances:
        if isinstance(instance, mol.ModificationInstance):
            mod_instances.add(instance)
        elif isinstance(instance, mol.AssociationInstance):
            assoc_instances.add(instance)
        elif isinstance(instance, mol.LocalizationInstance):
            if loc:
                raise AssertionError
            loc = instance

    return mol.MoleculeInstance(mol_def, mod_instances, assoc_instances, loc)


# CONVERTING CONTINGENCIES TO SETS OF STATES
def set_of_states_from_contingencies(contingencies: tg.List[con.Contingency]) -> venn.Set:
    if not contingencies:
        return venn.UniversalSet()

    for contingency in contingencies:
        assert contingency.target == contingencies[0].target
        assert contingency.type in [con.ContingencyType.inhibition, con.ContingencyType.requirement]

    requirements = []
    inhibitions = []

    for contingency in contingencies:
        if contingency.type == con.ContingencyType.requirement:
            requirements.append(set_of_states_from_effector(contingency.effector))

        elif contingency.type == con.ContingencyType.inhibition:
            inhibitions.append(set_of_states_from_effector(contingency.effector))

    required_set = venn.nested_expression_from_list_and_binary_op(requirements, venn.Intersection)
    inhibited_set = venn.Complement(venn.nested_expression_from_list_and_binary_op(inhibitions, venn.Union))

    if requirements and inhibitions:
        return venn.Intersection(required_set, inhibited_set)

    elif inhibitions:
        return inhibited_set

    elif requirements:
        return required_set


def source_set_of_states_from_reaction(reaction: rxn.Reaction) -> venn.Set:
    source_state = reaction.source
    product_state = reaction.product

    if not source_state and product_state:

        return venn.Complement(venn.PropertySet(product_state))

    elif source_state and not product_state:
        return venn.PropertySet(source_state)

    elif source_state and product_state:
        return venn.Intersection(venn.Complement(venn.PropertySet(product_state)),
                                 venn.PropertySet(source_state))

    else:
        raise AssertionError


def set_of_states_from_effector(effector: eff.Effector) -> venn.Set:
    if isinstance(effector, eff.StateEffector):
        return venn.PropertySet(effector.expr)

    elif isinstance(effector, eff.NotEffector):
        return venn.Complement(set_of_states_from_effector(effector.expr))

    elif isinstance(effector, eff.AndEffector):
        return venn.Intersection(set_of_states_from_effector(effector.left_expr), set_of_states_from_effector(effector.right_expr))

    elif isinstance(effector, eff.OrEffector):
        return venn.Union(set_of_states_from_effector(effector.left_expr), set_of_states_from_effector(effector.right_expr))

    else:
        raise AssertionError


# PROTECTED HELPERS
def _instances(mol_def: mol.MoleculeDefinition, state: sta.State, negate: bool) -> tg.List[mol.Instance]:
    if isinstance(state, sta.CovalentModificationState):
        matching_defs = [x for x in mol_def.modification_defs if x.spec.is_subspecification_of(state.substrate)]
        matching_instances = []
        for matching_def in matching_defs:
            if not negate:
                matching_instances.append(mol.ModificationInstance(matching_def,
                                                                   _molecule_modifier_from_state_modifier(state.modifier)))
            else:
                matching_instances.extend(mol.ModificationInstance(matching_def,
                                                                   _molecule_modifier_from_state_modifier(state.modifier)).complementary_instances())

        return matching_instances

    elif isinstance(state, sta.InterProteinInteractionState) or isinstance(state, sta.IntraProteinInteractionState):
        first_defs = [assoc_def for assoc_def in mol_def.association_defs
                      for valid_partner_spec in assoc_def.valid_partners
                      if assoc_def.spec.is_subspecification_of(state.first_component) and
                      valid_partner_spec.is_subspecification_of(state.second_component)]
        assert len(first_defs) <= 1
        matching_instances = []
        for matching_def in first_defs:
            partners = [spec for spec in matching_def.valid_partners if spec.is_subspecification_of(state.second_component)]
            assert len(partners) == 1

            if not negate:
                matching_instances.append(mol.AssociationInstance(matching_def, mol.OccupationStatus.occupied_known_partner,
                                                                  partners[0]))
            else:
                matching_instances.extend(mol.AssociationInstance(matching_def, mol.OccupationStatus.occupied_known_partner,
                                                                  partners[0]).complementary_instances())

        second_defs = [assoc_def for assoc_def in mol_def.association_defs
                       for valid_partner_spec in assoc_def.valid_partners
                       if assoc_def.spec.is_subspecification_of(state.second_component) and
                       valid_partner_spec.is_subspecification_of(state.first_component)]
        assert len(second_defs) <= 1
        for matching_def in second_defs:
            partners = [spec for spec in matching_def.valid_partners if spec.is_subspecification_of(state.first_component)]
            assert len(partners) == 1

            if not negate:
                matching_instances.append(mol.AssociationInstance(matching_def, mol.OccupationStatus.occupied_known_partner,
                                                                  partners[0]))
            else:
                matching_instances.extend(mol.AssociationInstance(matching_def, mol.OccupationStatus.occupied_known_partner,
                                                                  partners[0]).complementary_instances())

    else:
        raise NotImplementedError

    return matching_instances


def _mod_spec_domain_from_state(state: sta.CovalentModificationState, reaction: rxn.Reaction):

    spec = state.substrate
    if not spec.residue:
        spec.residue = _kinase_residue_name(reaction.subject)

    return state


def _mod_def_from_state_and_reaction(state: sta.CovalentModificationState, reaction: rxn.Reaction):

    state = _mod_spec_domain_from_state(state, reaction)
    mod_def = mol.ModificationDefinition(state.substrate,
                                         {mol.Modifier.unmodified, _molecule_modifier_from_state_modifier(state.modifier)})

    return mod_def


def _assoc_spec_domain_from_state(state: tg.Union[sta.InterProteinInteractionState, sta.IntraProteinInteractionState]):
    first_spec = state.first_component
    second_spec = state.second_component

    if not first_spec.domain:
        first_spec.domain = _assoc_domain_from_partner_spec(state.second_component)

    if not second_spec.domain:
        second_spec.domain = _assoc_domain_from_partner_spec(state.first_component)

    return first_spec, second_spec


def _assoc_defs_from_state(state: tg.Union[sta.InterProteinInteractionState, sta.IntraProteinInteractionState]):
    first_spec, second_spec = _assoc_spec_domain_from_state(state)

    first_def = mol.AssociationDefinition(first_spec, {second_spec})
    second_def = mol.AssociationDefinition(second_spec, {first_spec})

    return first_def, second_def


def _kinase_residue_name(spec: spe.Specification) -> str:
    return '{0}site'.format(spec.name)


def _assoc_domain_from_partner_spec(spec: spe.Specification) -> str:
    return '{0}assoc'.format(spec.name)


def _molecule_modifier_from_state_modifier(state_mod: sta.StateModifier) -> mol.Modifier:
    if state_mod == sta.StateModifier.phosphor:
        return mol.Modifier.phosphorylated

    elif state_mod == sta.StateModifier.ubiquitin:
        return mol.Modifier.ubiquitinated

    elif state_mod == sta.StateModifier.truncated:
        return mol.Modifier.trucated

    elif state_mod == sta.StateModifier.unmodified:
        raise NotImplementedError

    else:
        raise NotImplementedError


def _update_defs(defs: tg.Set[mol.Definition], new_def: mol.Definition):
    if isinstance(new_def, mol.AssociationDefinition):
        matches = lambda x, y: x == y
        update = lambda old_def, new_def: old_def.valid_partners.update(new_def.valid_partners)
    else:
        matches = lambda x, y: x.is_equivalent_to(y)
        update = lambda old_def, new_def: old_def.valid_modifiers.update(new_def.valid_modifiers)

    # if update is used the function will return a None which is stored in found_updateable_def resulting in an non-empty list
    found_updatable_def = [update(the_def, new_def) for the_def in defs if matches(the_def.spec, new_def.spec)]

    if not found_updatable_def:
        defs.add(new_def)


