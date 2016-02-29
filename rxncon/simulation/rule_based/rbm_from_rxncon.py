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






def mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(mol_def: mdf.MoleculeDefinition,
                                                                   reaction: rxn.Reaction,
                                                                   contingencies: tg.List[con.Contingency],
                                                                   disjunctify: bool=True) \
        -> tg.Tuple[mins.MoleculeInstance, mins.MoleculeInstance]:
    property_sets = mol_property_sets_from_mol_def_and_state_sets(mol_def, state_set_from_contingencies(contingencies), disjunctify)

    lhs_rhs_property_set_pairs = \
        [x for x in mol_property_pairs_from_mol_def_and_source_state_set(mol_def, source_state_set_from_reaction(reaction))
         if is_property_pair_valid_for_reaction(mol_def, x, reaction)]

    assert len(lhs_rhs_property_set_pairs) == 1
    lhs_source = lhs_rhs_property_set_pairs[0][0]
    rhs_source = lhs_rhs_property_set_pairs[0][1]

    pairs = []
    for property_set in property_sets:
        lhs_molecule = mfr.mol_instance_from_mol_def_and_property_set(mol_def, venn.Intersection(property_set, lhs_source))
        rhs_molecule = mfr.mol_instance_from_mol_def_and_property_set(mol_def, venn.Intersection(property_set, rhs_source))

        pairs.append((lhs_molecule, rhs_molecule))

    return pairs


def mol_property_sets_from_mol_def_and_state_sets(mol_def: mdf.MoleculeDefinition,
                                                  state_set: venn.Set,
                                                  disjunctify: bool=True) -> tg.List[venn.Set]:
    prop_sets = mfr.property_set_from_mol_def_and_state_set(mol_def, state_set).to_union_list_form()

    if disjunctify:
        return venn.gram_schmidt_disjunctify(prop_sets)

    else:
        return prop_sets


def mol_property_pairs_from_mol_def_and_source_state_set(mol_def: mdf.MoleculeDefinition, state_set: venn.Set):
    lhs_sets = mfr.property_set_from_mol_def_and_state_set(mol_def, state_set).to_union_list_form()
    rhs_sets = mfr.property_set_from_mol_def_and_state_set(mol_def, venn.Complement(state_set)).to_union_list_form()

    tuples = []
    for lhs in lhs_sets:
        for rhs in rhs_sets:
            tuples.append((lhs, rhs))

    return tuples


def is_property_pair_valid_for_reaction(mol_def: mdf.MoleculeDefinition,
                                        prop_tuple: tg.Tuple[venn.Set, venn.Set],
                                        reaction: rxn.Reaction) -> bool:
    lhs = prop_tuple[0].simplified_form()
    assert isinstance(lhs, venn.PropertySet)
    lhs_prop = lhs.value

    rhs = prop_tuple[1].simplified_form()
    assert isinstance(rhs, venn.PropertySet)
    rhs_prop = rhs.value

    return mfr.mol_def_and_property_match_state(mol_def, lhs_prop, reaction.source, negate=False) and\
        mfr.mol_def_and_property_match_state(mol_def, rhs_prop, reaction.product, negate=False)


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




