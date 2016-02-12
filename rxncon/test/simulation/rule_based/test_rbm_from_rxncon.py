import pytest

import rxncon.core.rxncon_system as rxs
import rxncon.semantics.molecule
import rxncon.semantics.molecule_from_rxncon
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.contingency as con
import rxncon.core.effector as eff

import rxncon.simulation.rule_based.rbm_from_rxncon as rfr
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.venntastic.sets as venn

import rxncon.input.quick.quick as qui

# @todo Fix the mutual exclusivity for associations and modifications.


# STRICT STATE SETS FROM STRICT CONTINGENCIES
def test_state_set_from_contingencies():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')

    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    b_ppi_e = rfs.reaction_from_string('B_ppi_E')
    b_dash_e = rfs.state_from_string('B--E')

    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    e_pplus = rfs.state_from_string("E-{P}")

    cont_b_dash_e = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(a_dash_b))  # B_ppi_E; ! A--B
    cont_e_pplus = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(e_pplus))  # B_ppi_E; ! E-{P}
    cont_a_ppi_b = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))  # A_ppi_B; x A--C

    rxncon = rxs.RxnConSystem([b_ppi_e, a_ppi_b, a_ppi_c, b_pplus_e], [cont_e_pplus, cont_b_dash_e, cont_a_ppi_b])

    strict_contingencies_state_set_b_ppi_e = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))

    strict_contingencies_state_set_a_ppi_b = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[1]))

    strict_contingencies_state_set_a_ppi_c = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[2]))

    strict_contingencies_state_set_b_pplus_e = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[3]))

    expected_b_ppi_e_strict_cont = venn.Intersection(venn.PropertySet(e_pplus), venn.PropertySet(a_dash_b))
    assert strict_contingencies_state_set_b_ppi_e.is_equivalent_to(expected_b_ppi_e_strict_cont)

    expected_a_ppi_b_strict_cont = venn.Intersection(venn.UniversalSet(),venn.Complement(venn.PropertySet(a_dash_c)))
    assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(expected_a_ppi_b_strict_cont)
    assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(venn.Complement(venn.PropertySet(a_dash_c)))

    expected_a_ppi_c_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_a_ppi_c.is_equivalent_to(expected_a_ppi_c_strict_cont)

    expected_b_pplus_e_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_b_pplus_e.is_equivalent_to(expected_b_pplus_e_strict_cont)


def test_state_set_from_contingencies_from_AND_complex():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")

    a_dash_c = rfs.state_from_string("A--C")
    c_dash_e = rfs.state_from_string("C--E")
    b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    strict_cont_state_set = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))


    assert strict_cont_state_set.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(a_dash_c), venn.PropertySet(c_dash_e)),
                                                                    venn.PropertySet(b_dash_f)))


# todo: this does not belong here ####################
def test_specification_set_from_state_set_single_ppi_no_contingency():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    rxncon = rxs.RxnConSystem([a_ppi_b], [])
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)

    strict_cont_state_set = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_spec_set = mol_defs["A"].specification_set_from_state_set(strict_cont_state_set)

    assert strict_spec_set.is_equivalent_to(venn.UniversalSet())


def test_specification_set_from_state_set_single_requirement_related():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(a_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont])
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)

    strict_cont_state_set = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_spec_set = mol_defs["A"].specification_set_from_state_set(strict_cont_state_set)

    assoc_def = [assoc_def for assoc_def in mol_defs["A"].association_defs if assoc_def.domain == "AssC"]

    assert strict_spec_set.is_equivalent_to(venn.PropertySet(
        rxncon.semantics.molecule.AssociationSpecification(assoc_def[0], rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)))


def test_specification_set_from_state_set_single_inhibition_related():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont])
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)

    strict_cont_state_set = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_spec_set = mol_defs["A"].specification_set_from_state_set(strict_cont_state_set)

    assoc_def = [assoc_def for assoc_def in mol_defs["A"].association_defs if assoc_def.domain == "AssC"]

    assert strict_spec_set.is_equivalent_to(venn.PropertySet(
        rxncon.semantics.molecule.AssociationSpecification(assoc_def[0], rxncon.semantics.molecule.OccupationStatus.not_occupied)))


def test_specification_set_from_state_set_single_requirement_not_related():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    b_ppi_c = rfs.reaction_from_string('B_ppi_C')
    b_dash_c = rfs.state_from_string('B--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(b_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, b_ppi_c], [cont])
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)

    strict_cont_state_set = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_spec_set = mol_defs["A"].specification_set_from_state_set(strict_cont_state_set)

    assert strict_spec_set.is_equivalent_to(venn.UniversalSet())


def test_specification_set_from_state_set_single_inhibited_not_related():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    b_ppi_c = rfs.reaction_from_string('B_ppi_C')
    b_dash_c = rfs.state_from_string('B--C')

    cont = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(b_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, b_ppi_c], [cont])
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)

    strict_cont_state_set = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_spec_set = mol_defs["A"].specification_set_from_state_set(strict_cont_state_set)

    assert strict_spec_set.is_equivalent_to(venn.UniversalSet())


def test_specification_set_from_state_set_quantitative_contingencies():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont_pos = con.Contingency(a_ppi_b, con.ContingencyType.positive, eff.StateEffector(a_dash_c))
    cont_neg = con.Contingency(a_ppi_b, con.ContingencyType.positive, eff.StateEffector(a_dash_c))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont_pos, cont_neg])
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)

    strict_cont_state_set = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_spec_set_A = mol_defs["A"].specification_set_from_state_set(strict_cont_state_set)
    strict_spec_set_C = mol_defs["C"].specification_set_from_state_set(strict_cont_state_set)

    assert strict_spec_set_A.is_equivalent_to(venn.UniversalSet())
    assert strict_spec_set_C.is_equivalent_to(venn.UniversalSet())


def test_specification_set_from_state_set_inhibition_one_molecule_not_related_to_both():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    cont_a_ppi_b = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))  # A_ppi_B; x A--C

    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c], [cont_a_ppi_b])
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)

    strict_cont_state_set_a_ppi_b = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))
    strict_spec_set_a_ppi_b = mol_defs["A"].specification_set_from_state_set(strict_cont_state_set_a_ppi_b)

    assoc_def = [assoc_def for assoc_def in mol_defs["A"].association_defs if assoc_def.domain == "AssC"]
    assert strict_spec_set_a_ppi_b.is_equivalent_to(venn.PropertySet(
        rxncon.semantics.molecule.AssociationSpecification(assoc_def[0], rxncon.semantics.molecule.OccupationStatus.not_occupied)))

    strict_spec_set_a_ppi_b = mol_defs["B"].specification_set_from_state_set(strict_cont_state_set_a_ppi_b)
    assert strict_spec_set_a_ppi_b.is_equivalent_to(venn.UniversalSet())


def test_specifictation_set_from_state_set_two_contingencies():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')

    b_ppi_e = rfs.reaction_from_string('B_ppi_E')

    b_pplus_e = rfs.reaction_from_string('E_p+_B')
    b_pplus = rfs.state_from_string("B-{P}")

    cont_b_dash_e = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(a_dash_b))  # B_ppi_E; ! A--B
    cont_e_pplus = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(b_pplus))  # B_ppi_E; ! B-{P}

    rxncon = rxs.RxnConSystem([b_ppi_e, a_ppi_b, b_pplus_e], [cont_e_pplus, cont_b_dash_e])

    strict_cont_state_set_b_ppi_e = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))

    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)

    strict_spec_set_b_ppi_e = mol_defs["B"].specification_set_from_state_set(strict_cont_state_set_b_ppi_e)

    assoc_def = [assoc_def for assoc_def in mol_defs["B"].association_defs if assoc_def.domain == "AssA"]

    assert strict_spec_set_b_ppi_e.is_equivalent_to(venn.Intersection(venn.PropertySet(
        rxncon.semantics.molecule.AssociationSpecification(assoc_def[0], rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)),
                                                                      venn.PropertySet(
                                                                          rxncon.semantics.molecule.ModificationSpecification(mol_defs["B"].modification_defs[0], "p"))))


def test_specifictation_set_from_state_boolean_complex():
    # todo throw something if state is not produced
    quick = qui.Quick("""
                        A_p+_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F
                        A_ppi_C
                        C_ppi_E
                        B_ppi_F""")
    rxncon = quick.rxncon_system
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(rxncon)
    strict_cont_state_set = rxncon.semantics.molecule_from_rxncon.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon.reactions[0]))

    strict_spec_set_B = mol_defs["B"].specification_set_from_state_set(strict_cont_state_set)
    strict_spec_set_A = mol_defs["A"].specification_set_from_state_set(strict_cont_state_set)
    strict_spec_set_C = mol_defs['C'].specification_set_from_state_set(strict_cont_state_set)

    assert strict_spec_set_B.is_equivalent_to(venn.PropertySet(
        rxncon.semantics.molecule.AssociationSpecification(mol_defs["B"].association_defs[0], rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)))
    assert strict_spec_set_A.is_equivalent_to(venn.PropertySet(
        rxncon.semantics.molecule.AssociationSpecification(mol_defs["A"].association_defs[0], rxncon.semantics.molecule.OccupationStatus.occupied_known_partner)))

    # @todo check that the association domains are the correct ones
    # @todo check the stuff for C and E
    print(strict_spec_set_C)

# todo: END####################


def test_simple_rxncon_system(simple_system):
    mol_defs = rxncon.semantics.molecule_from_rxncon.molecule_defs_from_rxncon(simple_system)

    mol_def_X = mol_defs['X']

    lhs_set = venn.Intersection(venn.PropertySet(rfs.state_from_string('X-{p}')),
                                venn.Complement(venn.PropertySet(rfs.state_from_string('X--Y'))))

    spec_set = mol_def_X.specification_set_from_state_set(lhs_set)

    spec_sets_overlapping = spec_set.to_union_list_form()
    spec_sets_disjunct = venn.gram_schmidt_disjunctify(spec_sets_overlapping)

    for ss in spec_sets_disjunct:
        print(mol_def_X.specification_from_specification_set(ss))


def test_single_rule():
    c_phos_a = rfs.reaction_from_string('C_p+_A_[x]')
    a_phos = rfs.state_from_string('A_[x]-{p}')
    a_phos_b = rfs.reaction_from_string('A_p+_B')

    cont = con.Contingency(a_phos_b, con.ContingencyType.requirement, eff.StateEffector(a_phos))

    rxncon = rxs.RxnConSystem([c_phos_a, a_phos_b], [cont])

    # print(rfr.molecule_defs_from_rxncon(rxncon)['B'].modification_defs[0].matching_state)

    rules = rfr.rules_from_rxncon(rxncon)

    for rule in rules:
        print()
        print(rule)


def test_single_ppi_rule():
    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    #c_ppi_f = rfs.reaction_from_string('C_ppi_F')
    b_ppi_e = rfs.reaction_from_string('B_ppi_E')
    a_dash_dash_c = rfs.state_from_string('A--C')
    b_dash_dash_e = rfs.state_from_string('B--E')
    #c_dash_dash_f = rfs.state_from_string('C--F')
    # A--B
    cont_A = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(a_dash_dash_c))
    cont_B = con.Contingency(a_ppi_b, con.ContingencyType.requirement, eff.StateEffector(b_dash_dash_e))
    #cont_C = con.Contingency(a_ppi_c, con.ContingencyType.requirement, eff.StateEffector(c_dash_dash_f))
    rxncon = rxs.RxnConSystem([a_ppi_b, a_ppi_c, b_ppi_e], [cont_A, cont_B])

    rules = rfr.rules_from_rxncon(rxncon)

@pytest.fixture
def simple_system():
    phosphorylation_reaction_1 = rfs.reaction_from_string('A_p+_X')
    phosphorylation_reaction_2 = rfs.reaction_from_string('B_p+_X')
    phosphorylation_reaction_3 = rfs.reaction_from_string('C_p+_X')

    binding_reaction = rfs.reaction_from_string('X_ppi_Y')
    phosphorylated_state = rfs.state_from_string('X-{p}')

    binding_contingency = con.Contingency(binding_reaction,
                                          con.ContingencyType.requirement,
                                          eff.StateEffector(phosphorylated_state))

    reactions = [phosphorylation_reaction_1, phosphorylation_reaction_2, phosphorylation_reaction_3, binding_reaction]
    contingencies = [binding_contingency]

    return rxs.RxnConSystem(reactions, contingencies)
