import typing as tg
import itertools as itt

import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.reaction as rxn
import rxncon.core.rxncon_system as rxs
import rxncon.venntastic.sets as venn
import rxncon.semantics.molecule_instance_from_rxncon as mfr
import rxncon.semantics.molecule_definition_from_rxncon as mdr
import rxncon.semantics.molecule_instance as mins
import rxncon.semantics.molecule_definition as mdf


class RuleBasedModelSupervisor:
    def __init__(self, rxncon: rxs.RxnConSystem):
        self.rxncon = rxncon
        self.mol_defs = mdr.MoleculeDefinitionSupervisor(rxncon).molecule_definitions

    def _generate_rule_prototypes(self):
        for reaction in self.rxncon.reactions:
            for mol_def in self.mol_defs.values():
                lhs_property_set = \
                    mfr.property_set_from_mol_def_and_state_set(mol_def, source_state_set_from_reaction(reaction))
                rhs_property_set = \
                    mfr.property_set_from_mol_def_and_state_set(mol_def, venn.Complement(source_state_set_from_reaction(reaction)))

                valid_lhs_to_rhs_pairs = _valid_property_set_pairs_for_mol_def_and_reaction(mol_def,
                                                                                            reaction,
                                                                                            lhs_property_set,
                                                                                            rhs_property_set)

                contingency_property_set = mfr.property_set_from_mol_def_and_state_set(
                    mol_def,
                    state_set_from_contingencies(self.rxncon.strict_contingencies_for_reaction(reaction))
                )





class RulePrototype:
    def __init__(self, lhs_molecules: tg.List[mins.MoleculeInstance], rhs_molecules: tg.List[mins.MoleculeInstance]):
        self.lhs_molecules = lhs_molecules
        self.rhs_molecules = rhs_molecules


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

    required_set = venn.nested_expression_from_list_and_binary_op(requirements, venn.Intersection)
    inhibited_set = venn.Complement(venn.nested_expression_from_list_and_binary_op(inhibitions, venn.Union))

    if requirements and inhibitions:
        return venn.Intersection(required_set, inhibited_set)

    elif inhibitions:
        return inhibited_set

    elif requirements:
        return required_set


def source_state_set_from_reaction(reaction: rxn.Reaction) -> venn.Set:
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


def _valid_property_set_pairs_for_mol_def_and_reaction(mol_def: mdf.MoleculeDefinition,
                                                       reaction: rxn.Reaction,
                                                       lhs: venn.Set,
                                                       rhs: venn.Set) -> tg.Tuple[venn.Set, venn.Set]:
    # These lists contain molecule property sets, we only want to keep the (lhs, rhs) pairs that are consistent with
    # the source/product states as given in the reaction definition.
    lhs_terms = lhs.to_union_list_form()
    rhs_terms = rhs.to_union_list_form()

    all_pairs = itt.product(lhs_terms, rhs_terms)

    valid_pairs = []

    for pair in all_pairs:
        if pair[0] and pair[1]:
            if mfr.mol_def_and_property_match_state(mol_def, pair[0], reaction.source, False) and \
               mfr.mol_def_and_property_match_state(mol_def, pair[1], reaction.product, False):
                valid_pairs.append(pair)
        elif pair[0] and not pair[1]:
            if mfr.mol_def_and_property_match_state(mol_def, pair[0], reaction.source, False) and not reaction.product:
                valid_pairs.append(pair)
        elif not pair[0] and pair[1]:
            if mfr.mol_def_and_property_match_state(mol_def, pair[1], reaction.product, False) and not reaction.source:
                valid_pairs.append(pair)

    return valid_pairs





