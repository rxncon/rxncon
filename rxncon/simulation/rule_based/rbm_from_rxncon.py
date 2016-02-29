import typing as tg
import itertools as itt

import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.core.reaction as rxn
import rxncon.core.rxncon_system as rxs
import rxncon.core.specification as spe
import rxncon.venntastic.sets as venn
import rxncon.semantics.molecule_instance_from_rxncon as mfr
import rxncon.semantics.molecule_definition_from_rxncon as mdr
import rxncon.semantics.molecule_instance as mins
import rxncon.semantics.molecule_definition as mdf
import rxncon.simulation.rule_based.rule_based_model as rbm


class RuleBasedModelSupervisor:
    def __init__(self, rxnconsys: rxs.RxnConSystem, disjunctify: bool=True):
        self.rxnconsys = rxnconsys
        self.mol_defs = mdr.MoleculeDefinitionSupervisor(self.rxnconsys).molecule_definitions
        self.rules = []  # type: tg.List[rbm.Rule]
        self.rule_based_model = None
        self._generate_rules(disjunctify)
        self._construct_rule_based_model()

        self._validate()

    def _generate_rules(self, disjunctify: bool):
        for reaction in self.rxnconsys.reactions:
            strict_conts = self.rxnconsys.strict_contingencies_for_reaction(reaction)
            involved_molecules = involved_molecule_specs_for_reaction_and_contingencies(reaction, strict_conts)
            mol_instance_pairs = []

            for mol_def in self.mol_defs:
                if mol_def.spec in involved_molecules:
                    mol_instance_pairs.append(mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(mol_def, reaction, strict_conts, disjunctify))

            rules_of_mol_instances = lhs_rhs_product(mol_instance_pairs)
            for rule_of_mol_instance in rules_of_mol_instances:
                lhs_reactants = reactants_from_molecule_instances(rule_of_mol_instance[0])
                rhs_reactants = reactants_from_molecule_instances(rule_of_mol_instance[1])
                self.rules.append(rbm.Rule(lhs_reactants, rhs_reactants, arrow_from_reaction(reaction), None))

    def _construct_rule_based_model(self):
        self.rule_based_model = rbm.RuleBasedModel(self.mol_defs, self.rules, None, None)

    def _validate(self):
        pass


def involved_molecule_specs_for_reaction_and_contingencies(reaction: rxn.Reaction,
                                                           strict_cont: tg.List[con.Contingency]) -> tg.List[spe.Specification]:
    pass


def reactants_from_molecule_instances(molecules: tg.List[mins.MoleculeInstance]) -> tg.List[rbm.Reactant]:
    pass


def lhs_rhs_product(reaction_molecules: tg.List[tg.List[tg.Tuple[mins.MoleculeInstance, mins.MoleculeInstance]]])\
        -> tg.List[tg.Tuple[tg.List[mins.MoleculeInstance], tg.List[mins.MoleculeInstance]]]:
    mol_product = itt.product(reaction_molecules)

    lhs_rhs = []

    for x in mol_product:
        lhs_rhs.append(tuple(zip(*x)))

    return lhs_rhs


def mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(mol_def: mdf.MoleculeDefinition,
                                                                   reaction: rxn.Reaction,
                                                                   contingencies: tg.List[con.Contingency],
                                                                   disjunctify: bool=True) \
        -> tg.List[tg.Tuple[mins.MoleculeInstance, mins.MoleculeInstance]]:
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


def mol_property_pairs_from_mol_def_and_source_state_set(mol_def: mdf.MoleculeDefinition, state_set: venn.Set) \
        -> tg.Tuple[venn.PropertySet, venn.PropertySet]:
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


