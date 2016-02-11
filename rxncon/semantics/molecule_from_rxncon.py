import functools as ft
import itertools as itt
import typing as tg
from collections import defaultdict

import rxncon.semantics.molecule as mol
import rxncon.core.rxncon_system as rxs
import rxncon.core.state as sta
import rxncon.core.reaction as rxn
import rxncon.core.specification as spe


class MoleculeDefinitionSupervisor:
    def __init__(self, rxnconsys: rxs.RxnConSystem):
        self.rxnconsys = rxnconsys
        self.molecule_definitions = {}
        self._generate_molecule_definitions()

    def molecule_definition_for_name(self, name: str):
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

                    assoc_defs = _state_to_assoc_defs(state)
                    _update_defs(name_to_assoc_defs[state.first_component.name], assoc_defs[0])
                    _update_defs(name_to_assoc_defs[state.second_component.name], assoc_defs[1])

                elif isinstance(state, sta.TranslocationState):
                    names.add(state.substrate.name)
                    names.add(reaction.subject.name)
                    name_to_locs[state.substrate.name].add(state.compartment)

                else:
                    raise NotImplementedError

        for name in names:
            self.molecule_definitions[name] = mol.MoleculeDefinition(name,
                                                                     name_to_mod_defs[name],
                                                                     name_to_assoc_defs[name],
                                                                     mol.LocalizationDefinition(name_to_locs[name]))


    def specification_set_from_state_set(self, state_set: venn.Set) -> venn.Set:
        # molecule_def + set of states --> set of specifications
        assert len(state_set.to_nested_list_form()) == 1
        spec_sets = []

        if state_set == venn.UniversalSet():
            return venn.UniversalSet()

        for single_property in state_set.to_nested_list_form()[0]:
            if isinstance(single_property, venn.PropertySet):
                state = single_property.value
                negate = False
            elif isinstance(single_property, venn.Complement):
                assert isinstance(single_property.expr, venn.PropertySet)
                state = single_property.expr.value
                negate = True
            else:
                raise AssertionError

            matching_specs = []

            for definition in self.modification_defs + self.association_defs:
                matching_specs += definition.match_with_state(state, negate=negate)

            matching_specs += self.localization_def.match_with_state(state, negate=negate)

            if matching_specs:
                spec_sets.append(venn.nested_expression_from_list_and_binary_op([venn.PropertySet(x) for x in matching_specs], venn.Union))
            else:
                spec_sets.append(venn.UniversalSet())


        if all(x == venn.EmptySet() for x in spec_sets):
            return venn.UniversalSet()

        return venn.nested_expression_from_list_and_binary_op(spec_sets, venn.Intersection)

    def specification_from_specification_set(self, spec_set: venn.Set) -> 'MoleculeInstance':
        spec_lists = remove_complements_from_spec_set(spec_set).to_nested_list_form()

        assert len(spec_lists) == 1
        spec_list = spec_lists[0]

        assert all(isinstance(x, venn.PropertySet) for x in spec_list)
        ass_specs = [spec.value for spec in spec_list if isinstance(spec.value, AssociationInstance)]
        mod_specs = [spec.value for spec in spec_list if isinstance(spec.value, ModificationInstance)]
        loc_specs = [spec.value for spec in spec_list if isinstance(spec.value, LocalizationInstance)]

        if len(loc_specs) == 1:
            loc_spec = loc_specs[0]
        elif len(loc_specs) == 0:
            loc_spec = None
        else:
            raise AssertionError

        return MoleculeInstance(self, mod_specs, ass_specs, loc_spec)



def state_set_from_contingencies(contingencies: tg.List[con.Contingency]) -> venn.Set:
    if not contingencies:
        return venn.UniversalSet()

    for contingency in contingencies:
        assert contingency.target == contingencies[0].target
        assert contingency.type in [con.ContingencyType.inhibition, con.ContingencyType.requirement]

    requirements = []
    inhibitions = []

    for contingency in contingencies:
        if contingency.type == con.ContingencyType.requirement:
            requirements.append(state_set_from_effector(contingency.effector))

        elif contingency.type == con.ContingencyType.inhibition:
            inhibitions.append(state_set_from_effector(contingency.effector))

    required = venn.nested_expression_from_list_and_binary_op(requirements, venn.Intersection)
    inhibited = venn.nested_expression_from_list_and_binary_op(inhibitions, venn.Union)

    if required != venn.EmptySet() and inhibited != venn.EmptySet():
        return venn.Intersection(required, venn.Complement(inhibited))

    elif required == venn.EmptySet():
        return venn.Complement(inhibited)

    elif inhibited == venn.EmptySet():
        return required


def contingency_configurations_from_quantitative_contingencies(contingencies: tg.List[con.Contingency]) -> tg.List[con.Contingency]:
    if not contingencies:
        return []

    reaction = contingencies[0].target

    for contingency in contingencies:
        assert contingency.target == reaction
        assert contingency.type in [con.ContingencyType.positive, con.ContingencyType.negative]

    config_pairs = [(effector, eff.NotEffector(effector)) for effector in [contingency.effector for contingency in contingencies]]

    configs = itt.product(*config_pairs)

    contingency_configs = []

    for config in configs:
        if len(config) == 1:
            total_effector = config[0]
        else:
            total_effector = ft.reduce(eff.AndEffector, config[1:], config[0])

        contingency_configs.append(con.Contingency(reaction, con.ContingencyType.requirement, total_effector))

    return contingency_configs


def source_state_set_from_reaction(reaction: rxn.Reaction) -> venn.Set:
    source_state = reaction.source
    product_state = reaction.product

    if not source_state and product_state:
        return venn.Complement(venn.PropertySet(product_state))

    elif source_state and not product_state:
        return venn.PropertySet(source_state)

    elif source_state and product_state:
        return venn.Intersection(venn.Complement(venn.PropertySet(product_state)), venn.PropertySet(source_state))

    else:
        raise AssertionError


def state_set_from_effector(effector: eff.Effector) -> venn.Set:
    if isinstance(effector, eff.StateEffector):
        return venn.PropertySet(effector.expr)

    elif isinstance(effector, eff.NotEffector):
        return venn.Complement(state_set_from_effector(effector.expr))

    elif isinstance(effector, eff.AndEffector):
        return venn.Intersection(state_set_from_effector(effector.left_expr), state_set_from_effector(effector.right_expr))

    elif isinstance(effector, eff.OrEffector):
        return venn.Union(state_set_from_effector(effector.left_expr), state_set_from_effector(effector.right_expr))

    else:
        raise AssertionError


def remove_complements_from_spec_set(spec_set):
    if isinstance(spec_set, venn.Intersection):
        return venn.Intersection(remove_complements_from_spec_set(spec_set.left_expr), remove_complements_from_spec_set(spec_set.right_expr))

    elif isinstance(spec_set, venn.Union):
        return venn.Union(remove_complements_from_spec_set(spec_set.left_expr), remove_complements_from_spec_set(spec_set.right_expr))

    elif isinstance(spec_set, venn.PropertySet):
        return spec_set

    elif isinstance(spec_set, venn.Complement):
        assert isinstance(spec_set.expr, venn.PropertySet)

        if not spec_set.expr.value:
            return venn.UniversalSet()

        complement_terms = [venn.PropertySet(x) for x in spec_set.expr.value.complementary_instances()]

        return venn.nested_expression_from_list_and_binary_op(complement_terms, venn.Union)

    else:
        raise AssertionError


def molecule_defs_from_rxncon(rxnconsys: rxs.RxnConSystem) -> tg.Dict[str, rxncon.semantics.molecule.MoleculeDefinition]:


def _mod_def_from_state_and_reaction(state: sta.CovalentModificationState, reaction: rxn.Reaction):
    spec = state.substrate

    if not spec.residue:
        spec.residue = _kinase_residue_name(reaction.subject)

    mod_def = mol.ModificationDefinition(spec,
                                         {mol.Modifier.unmodified, _molecule_modifier_from_state_modifier(state.modifier)})

    return mod_def


def _state_to_assoc_defs(state: tg.Union[sta.InterProteinInteractionState, sta.IntraProteinInteractionState]):
    first_spec = state.first_component
    second_spec = state.second_component

    if not first_spec.domain:
        first_spec.domain = _assoc_domain_from_partner_spec(state.second_component)

    if not second_spec.domain:
        second_spec.domain = _assoc_domain_from_partner_spec(state.first_component)

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
        update = lambda old_def, new_def: old_def.valid_partners.add(new_def.valid_partners)
    else:
        matches = lambda x, y: x.is_equivalent_to(y)
        update = lambda old_def, new_def: old_def.valid_modifiers.add(new_def.valid_modifiers)

    found_updatable_def = False
    for the_def in defs:
        if matches(the_def.spec, new_def.spec):
            update(the_def, new_def)
            found_updatable_def = True

    if not found_updatable_def:
        defs.add(new_def)








