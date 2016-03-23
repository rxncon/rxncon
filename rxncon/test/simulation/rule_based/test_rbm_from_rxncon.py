import pytest
from rxncon.core.rxncon_system import RxnConSystem
from rxncon.syntax.rxncon_from_string import reaction_from_string, state_from_string

import rxncon.core.specification as spec
import rxncon.core.effector as eff
import rxncon.simulation.rule_based.rbm_from_rxncon as rfr
from rxncon.core import contingency as con
from rxncon.input.quick import quick as qui
from rxncon.venntastic import sets as venn
import rxncon.semantics.molecule_definition_from_rxncon as mdr
import rxncon.semantics.molecule_instance as moi
import rxncon.simulation.rule_based.rule_based_model as rbm


# TEST RULE GENERATION FOR A LIST OF RXNCON SYSTEMS
def test_generate_rules(rxn_systems, expected_rules):
    for i, system in enumerate(rxn_systems):
        actual_rules = rfr.RuleBasedModelSupervisor(system).rules

        for rule in actual_rules:
            assert rule in expected_rules[i]

        for rule in expected_rules[i]:
            assert rule in actual_rules


# # MASTERTEST testing the lhs/rhs MoleculeInstances that appear.
# def test_mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(simple_system: RxnConSystem):
#     phosphorylation_reaction_X = simple_system.reactions[0]
#
#     mol_defs = mdr.MoleculeDefinitionSupervisor(simple_system).molecule_definitions
#     mol_def_X = mol_defs['X']
#
#     strict_conts = simple_system.strict_contingencies_for_reaction(phosphorylation_reaction_X)
#
#     phosphorylation_reaction_X_pairs = rfr.mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(mol_def_X, phosphorylation_reaction_X, strict_conts)
#     mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "Asite"]
#
#     expected_lhs = moi.MoleculeInstance(mol_def_X, {moi.ModificationPropertyInstance(mod_defs[0], moi.Modifier("u"))}, set(), None)
#     expected_rhs = moi.MoleculeInstance(mol_def_X, {moi.ModificationPropertyInstance(mod_defs[0], moi.Modifier("p"))}, set(), None)
#     assert phosphorylation_reaction_X_pairs == [(expected_lhs, expected_rhs)]
#
#     phosphorylation_reaction_X_at_resi = simple_system.reactions[1]
#     mol_defs = mdr.MoleculeDefinitionSupervisor(simple_system).molecule_definitions
#     mol_def_X = mol_defs['X']
#     strict_conts = simple_system.strict_contingencies_for_reaction(phosphorylation_reaction_X_at_resi)
#
#     phosphorylation_reaction_X_at_resi_pairs = rfr.mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(mol_def_X, phosphorylation_reaction_X_at_resi, strict_conts)
#
#     mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "r"]
#
#     expected_lhs_phosphorylation_reaction_X_at_resi = moi.MoleculeInstance(mol_def_X, {moi.ModificationPropertyInstance(mod_defs[0], moi.Modifier("u"))}, set(), None)
#     expected_rhs_phosphorylation_reaction_X_at_resi = moi.MoleculeInstance(mol_def_X, {moi.ModificationPropertyInstance(mod_defs[0], moi.Modifier("p"))}, set(), None)
#     assert phosphorylation_reaction_X_at_resi_pairs == [(expected_lhs_phosphorylation_reaction_X_at_resi,
#                                                          expected_rhs_phosphorylation_reaction_X_at_resi)]
#
#
# def test_set_of_states_from_single_state_effector(simple_system: RxnConSystem):
#     cont1 = simple_system.contingencies[0]
#     cont2 = simple_system.contingencies[1]
#
#     state_set_of_effector_cont1 = rfr.state_set_from_effector(cont1.effector)
#     state_set_of_effector_cont2 = rfr.state_set_from_effector(cont2.effector)
#
#     assert cont1.target == rxncon_from_string.reaction_from_string('Y_ppi_X_[d]')
#     assert cont1.type == con.ContingencyType.requirement
#     assert state_set_of_effector_cont1.is_equivalent_to(venn.PropertySet(rxncon_from_string.state_from_string('X-{P}')))
#
#     assert cont2.target == rxncon_from_string.reaction_from_string('Y_ppi_X_[d]')
#     assert cont2.type == con.ContingencyType.inhibition
#     assert state_set_of_effector_cont2.is_equivalent_to(venn.PropertySet(rxncon_from_string.state_from_string('A--X')))
#
#
# def test_set_of_states_from_nested_AND_effectors():
#     quick = qui.Quick("""A_ppi_B; ! <comp>
#                         <comp>; AND A--C
#                         <comp>; AND C--E
#                         <comp>; AND B--F""")
#
#     expected_a_dash_c = rxncon_from_string.state_from_string("A--C")
#     expected_c_dash_e = rxncon_from_string.state_from_string("C--E")
#     expected_b_dash_f = rxncon_from_string.state_from_string("B--F")
#
#     rxncon = quick.rxncon_system
#
#     set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
#
#     assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c),
#                                                                                           venn.PropertySet(expected_c_dash_e)),
#                                                                         venn.PropertySet(expected_b_dash_f)))
#
#
# def test_set_of_states_from_nested_OR_effectors():
#     quick = qui.Quick("""A_ppi_B; ! <comp>
#                         <comp>; OR A--C
#                         <comp>; OR C--E
#                         <comp>; OR B--F""")
#
#     expected_a_dash_c = rxncon_from_string.state_from_string("A--C")
#     expected_c_dash_e = rxncon_from_string.state_from_string("C--E")
#     expected_b_dash_f = rxncon_from_string.state_from_string("B--F")
#
#     rxncon = quick.rxncon_system
#
#     set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
#     assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Union(venn.PropertySet(expected_a_dash_c),
#                                                                                           venn.PropertySet(expected_c_dash_e)),
#                                                                         venn.PropertySet(expected_b_dash_f)))
#
#
# def test_set_of_states_from_nested_AND_OR_effector():
#     quick = qui.Quick("""A_ppi_B; ! <comp>
#                         <comp>; AND <c1>
#                         <comp>; AND <c2>
#                         <c1>; OR A--C
#                         <c1>; OR C--E
#                         <c2>; OR B--F
#                         <c2>; OR B--D""")
#
#     expected_a_dash_c = rxncon_from_string.state_from_string("A--C")
#     expected_c_dash_e = rxncon_from_string.state_from_string("C--E")
#     expected_b_dash_f = rxncon_from_string.state_from_string("B--F")
#     expected_b_dash_d = rxncon_from_string.state_from_string("B--D")
#
#     rxncon = quick.rxncon_system
#
#     set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
#     assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Union(venn.PropertySet(expected_a_dash_c),
#                                                                                           venn.PropertySet(expected_c_dash_e)),
#                                                                         venn.Union(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))
#
#
# def test_set_of_states_from_nested_OR_AND_effector():
#     quick = qui.Quick("""A_ppi_B; ! <comp>
#                         <comp>; OR <c1>
#                         <comp>; OR <c2>
#                         <c1>; AND A--C
#                         <c1>; AND C--E
#                         <c2>; AND B--F
#                         <c2>; AND B--D""")
#
#     expected_a_dash_c = rxncon_from_string.state_from_string("A--C")
#     expected_c_dash_e = rxncon_from_string.state_from_string("C--E")
#     expected_b_dash_f = rxncon_from_string.state_from_string("B--F")
#     expected_b_dash_d = rxncon_from_string.state_from_string("B--D")
#
#
#     rxncon = quick.rxncon_system
#
#     set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
#     assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
#                                                                                           venn.PropertySet(expected_c_dash_e)),
#                                                                         venn.Intersection(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))
#
#
# def test_set_of_states_from_complex_effector():
#     quick = qui.Quick("""A_ppi_B; ! <comp>
#                         <comp>; OR <c1>
#                         <comp>; OR <c2>
#                         <c1>; AND A--C
#                         <c1>; AND C--E
#                         <c2>; NOT <c3>
#                         <c3>; AND B--F
#                         <c3>; AND B--D""")
#
#     expected_a_dash_c = rxncon_from_string.state_from_string("A--C")
#     expected_c_dash_e = rxncon_from_string.state_from_string("C--E")
#     expected_b_dash_f = rxncon_from_string.state_from_string("B--F")
#     expected_b_dash_d = rxncon_from_string.state_from_string("B--D")
#
#     rxncon = quick.rxncon_system
#
#     set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
#     assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
#                                                                                           venn.PropertySet(expected_c_dash_e)),
#                                                                         venn.Complement(venn.Intersection(venn.PropertySet(expected_b_dash_d),
#                                                                                                           venn.PropertySet(expected_b_dash_f)))))
#
#
# def test_state_set_from_contingencies(simple_system):
#
#
#     a_ppi_b = rxncon_from_string.reaction_from_string('A_ppi_B')
#     a_dash_b = rxncon_from_string.state_from_string('A--B')
#
#     a_ppi_c = rxncon_from_string.reaction_from_string('A_ppi_C')
#     a_dash_c = rxncon_from_string.state_from_string('A--C')
#
#     b_ppi_e = rxncon_from_string.reaction_from_string('B_ppi_E')
#
#     b_pplus_e = rxncon_from_string.reaction_from_string('B_p+_E')
#     e_pplus = rxncon_from_string.state_from_string("E-{P}")
#
#     cont_b_dash_e = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(a_dash_b))  # B_ppi_E; ! A--B
#     cont_e_pplus = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(e_pplus))  # B_ppi_E; ! E-{P}
#     cont_a_ppi_b = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))  # A_ppi_B; x A--C
#
#     rxncon = RxnConSystem([b_ppi_e, a_ppi_b, a_ppi_c, b_pplus_e], [cont_e_pplus, cont_b_dash_e, cont_a_ppi_b])
#
#     strict_contingencies_state_set_b_ppi_e = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(b_ppi_e))
#     strict_contingencies_state_set_a_ppi_b = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_b))
#     strict_contingencies_state_set_a_ppi_c = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_c))
#     strict_contingencies_state_set_b_pplus_e = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(b_pplus_e))
#
#     expected_b_ppi_e_strict_cont = venn.Intersection(venn.PropertySet(e_pplus), venn.PropertySet(a_dash_b))
#     assert strict_contingencies_state_set_b_ppi_e.is_equivalent_to(expected_b_ppi_e_strict_cont)
#
#     expected_a_ppi_b_strict_cont = venn.Complement(venn.PropertySet(a_dash_c))
#     assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(expected_a_ppi_b_strict_cont)
#
#     expected_a_ppi_c_strict_cont = venn.UniversalSet()
#     assert strict_contingencies_state_set_a_ppi_c.is_equivalent_to(expected_a_ppi_c_strict_cont)
#
#     expected_b_pplus_e_strict_cont = venn.UniversalSet()
#     assert strict_contingencies_state_set_b_pplus_e.is_equivalent_to(expected_b_pplus_e_strict_cont)
#
#
# def test_set_of_states_from_contingencies_FOR_quant():
#     # todo
#     pass
#
#
# def test_set_of_states_from_strict_contingencies_AND():
#     quick = qui.Quick("""A_ppi_B; ! <comp>
#                         <comp>; AND A--C
#                         <comp>; AND C--E
#                         <comp>; AND B--F""")
#
#     expected_a_dash_c = rxncon_from_string.state_from_string("A--C")
#     expected_c_dash_e = rxncon_from_string.state_from_string("C--E")
#     expected_b_dash_f = rxncon_from_string.state_from_string("B--F")
#
#     rxncon = quick.rxncon_system
#
#     strict_cont_state_set = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rxncon_from_string.reaction_from_string('A_ppi_B')))
#
#     assert strict_cont_state_set.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c), venn.PropertySet(expected_c_dash_e)),
#                                                                     venn.PropertySet(expected_b_dash_f)))
#
#
# def test_source_state_set_from_reaction(simple_system: RxnConSystem):
#
#     phosphorylation_reaction_X = simple_system.reactions[0]
#     phosphorylation_reaction_X_at_resi = simple_system.reactions[1]
#     ubiquitination_reaction_X_at_resi = simple_system.reactions[2]
#     dephosphorylation_reaction_X = simple_system.reactions[3]
#     phosphortransfer_reaction_X = simple_system.reactions[4]
#     binding_reaction_A_X = simple_system.reactions[5]
#     binding_reaction_Y_X = simple_system.reactions[6]
#
#     assert rfr.source_state_set_from_reaction(phosphorylation_reaction_X).is_equivalent_to(venn.Complement(
#                                                                                                 venn.PropertySet(
#                                                                                                     rxncon_from_string.state_from_string(
#                                                                                                         'X-{P}'))))
#     assert rfr.source_state_set_from_reaction(phosphorylation_reaction_X_at_resi).is_equivalent_to(venn.Complement(
#                                                                                                         venn.PropertySet(
#                                                                                                             rxncon_from_string.state_from_string('X_(r)-{P}'))))
#     assert rfr.source_state_set_from_reaction(ubiquitination_reaction_X_at_resi).is_equivalent_to(venn.Complement(
#                                                                                                         venn.PropertySet(
#                                                                                                             rxncon_from_string.state_from_string('X_(r)-{Ub}'))))
#     assert rfr.source_state_set_from_reaction(dephosphorylation_reaction_X).is_equivalent_to(venn.PropertySet(rxncon_from_string.state_from_string('X-{P}')))
#
#     # todo: B_pt_E are two reactions in one E_p+_X -> X_[Eside] and X_p-_E -> E_[Xside]
#     # todo: B_[n]_apt_B_[m] auto phosphortransfer B is the same molecule B_[n]_p+_B_[m] -> B_[m] and B_[m]_p-_B_[n] -> B_B[n]
#     # todo: B_apt_B auto phosphortransfer B is the same molecule B_p+_B -> B_[Bsite1] and B_p-_B -> B_B[Site2]
#
#     assert rfr.source_state_set_from_reaction(phosphortransfer_reaction_X).is_equivalent_to(venn.Intersection(venn.Complement(venn.PropertySet(rxncon_from_string.state_from_string('X-{p}'))),
#                                                                                                               venn.PropertySet(rxncon_from_string.state_from_string('E-{p}'))))
#     assert rfr.source_state_set_from_reaction(binding_reaction_A_X).is_equivalent_to(venn.Complement(venn.PropertySet(rxncon_from_string.state_from_string('A--X_[d]'))))
#     assert rfr.source_state_set_from_reaction(binding_reaction_Y_X).is_equivalent_to(venn.Complement(venn.PropertySet(rxncon_from_string.state_from_string('Y--X_[d]'))))
#
#
# def test_mol_property_pairs_from_mol_def_and_source_state_set(simple_system):
#     mol_defs = mdr.MoleculeDefinitionSupervisor(simple_system).molecule_definitions
#     mol_def_X = mol_defs['X']
#
#     phosphorylation_reaction_X = simple_system.reactions[0]
#     phosphorylation_reaction_X_at_resi = simple_system.reactions[1]
#     ubiquitination_reaction_X_at_resi = simple_system.reactions[2]
#     dephosphorylation_reaction_X = simple_system.reactions[3]
#     phosphortransfer_reaction_X = simple_system.reactions[4]
#     binding_reaction_A_X = simple_system.reactions[5]
#     binding_reaction_Y_X = simple_system.reactions[6]
#
#
#     mol_prop_pair_phosphorylation_reaction_X = rfr.mol_property_pairs_from_mol_def_and_reaction(mol_def_X,
#                                                                                                 rfr.source_state_set_from_reaction(phosphorylation_reaction_X))
#     assert len(mol_prop_pair_phosphorylation_reaction_X) == 1
#     assert len(mol_prop_pair_phosphorylation_reaction_X[0]) == 2
#
#     mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "Asite"]
#     assert len(mod_defs) == 1
#     mod_def_Asite = mod_defs[0]
#     assert mol_prop_pair_phosphorylation_reaction_X[0][0] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_Asite, moi.Modifier.unmodified))
#     assert mol_prop_pair_phosphorylation_reaction_X[0][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_Asite, moi.Modifier.phosphorylated))
#
#
#     mol_prop_pair_phosphorylation_reaction_X_at_resi = rfr.mol_property_pairs_from_mol_def_and_reaction(mol_def_X,
#                                                                                                         rfr.source_state_set_from_reaction(phosphorylation_reaction_X_at_resi))
#
#     assert len(mol_prop_pair_phosphorylation_reaction_X_at_resi) == 2
#     assert len(mol_prop_pair_phosphorylation_reaction_X_at_resi[0]) == 2
#     assert len(mol_prop_pair_phosphorylation_reaction_X_at_resi[1]) == 2
#
#     mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "r"]
#     assert len(mod_defs) == 1
#     mod_def_r = mod_defs[0]
#
#     assert mol_prop_pair_phosphorylation_reaction_X_at_resi[0][0] != mol_prop_pair_phosphorylation_reaction_X_at_resi[1][0]
#     assert mol_prop_pair_phosphorylation_reaction_X_at_resi[0][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_r, moi.Modifier.phosphorylated))
#     assert mol_prop_pair_phosphorylation_reaction_X_at_resi[1][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_r, moi.Modifier.phosphorylated))
#
#
#     mol_prop_pair_ubiquitination_reaction_X_at_resi = rfr.mol_property_pairs_from_mol_def_and_reaction(mol_def_X,
#                                                                                                        rfr.source_state_set_from_reaction(ubiquitination_reaction_X_at_resi))
#
#     assert len(mol_prop_pair_ubiquitination_reaction_X_at_resi) == 2
#     assert len(mol_prop_pair_ubiquitination_reaction_X_at_resi[0]) == 2
#     assert len(mol_prop_pair_ubiquitination_reaction_X_at_resi[1]) == 2
#
#
#     assert mol_prop_pair_ubiquitination_reaction_X_at_resi[0][0] != mol_prop_pair_ubiquitination_reaction_X_at_resi[1][0]
#     assert mol_prop_pair_ubiquitination_reaction_X_at_resi[0][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_r, moi.Modifier.ubiquitinated))
#     assert mol_prop_pair_ubiquitination_reaction_X_at_resi[1][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_r, moi.Modifier.ubiquitinated))
#
#
#     mol_prop_pair_dephosphorylation_reaction_X = rfr.mol_property_pairs_from_mol_def_and_reaction(mol_def_X,
#                                                                                                   rfr.source_state_set_from_reaction(dephosphorylation_reaction_X))
#
#
#     assert len(mol_prop_pair_dephosphorylation_reaction_X) == 1
#     assert len(mol_prop_pair_dephosphorylation_reaction_X[0]) == 2
#
#     mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "Bsite"]
#     assert len(mod_defs) == 1
#     mod_def_Bsite = mod_defs[0]
#
#     assert mol_prop_pair_dephosphorylation_reaction_X[0][0] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_Bsite, moi.Modifier.phosphorylated))
#     assert mol_prop_pair_dephosphorylation_reaction_X[0][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_Bsite, moi.Modifier.unmodified))
#
#
#     mol_prop_pair_binding_reaction_A_X = rfr.mol_property_pairs_from_mol_def_and_reaction(mol_def_X,
#                                                                                           rfr.source_state_set_from_reaction(binding_reaction_A_X))
#
#     assert len(mol_prop_pair_binding_reaction_A_X) == 2
#     assert len(mol_prop_pair_binding_reaction_A_X[0]) == 2
#     assert len(mol_prop_pair_binding_reaction_A_X[1]) == 2
#
#     assoc_defs = [x for x in mol_def_X.association_defs if x.spec.domain == "d"]
#     assert len(assoc_defs) == 1
#     assoc_def_d = assoc_defs[0]
#
#     assert mol_prop_pair_binding_reaction_A_X[0][0] != mol_prop_pair_binding_reaction_A_X[1][0]
#     assert mol_prop_pair_binding_reaction_A_X[0][1] == venn.PropertySet(moi.AssociationPropertyInstance(assoc_def_d,
#                                                                                                         moi.OccupationStatus.occupied_known_partner,
#                                                                                                         spec.Specification('A', 'Xassoc', None, None)))
#     assert mol_prop_pair_binding_reaction_A_X[1][1] == venn.PropertySet(moi.AssociationPropertyInstance(assoc_def_d,
#                                                                                                         moi.OccupationStatus.occupied_known_partner,
#                                                                                                         spec.Specification('A', 'Xassoc', None, None)))
#     mol_prop_pair_binding_reaction_Y_X = rfr.mol_property_pairs_from_mol_def_and_reaction(mol_def_X,
#                                                                                           rfr.source_state_set_from_reaction(binding_reaction_Y_X))
#
#     assert len(mol_prop_pair_binding_reaction_Y_X) == 2
#     assert len(mol_prop_pair_binding_reaction_Y_X[0]) == 2
#     assert len(mol_prop_pair_binding_reaction_Y_X[1]) == 2
#
#     assert mol_prop_pair_binding_reaction_Y_X[0][0] != mol_prop_pair_binding_reaction_Y_X[1][0]
#     assert mol_prop_pair_binding_reaction_Y_X[0][1] == venn.PropertySet(moi.AssociationPropertyInstance(assoc_def_d,
#                                                                                                         moi.OccupationStatus.occupied_known_partner,
#                                                                                                         spec.Specification('Y', 'Xassoc', None, None)))
#     assert mol_prop_pair_binding_reaction_Y_X[1][1] == venn.PropertySet(moi.AssociationPropertyInstance(assoc_def_d,
#                                                                                                         moi.OccupationStatus.occupied_known_partner,
#                                                                                                         spec.Specification('Y', 'Xassoc', None, None)))
#     #is_property_pair_valid_for_reaction(mol_def, x, reaction)
#
#     mol_prop_pair_phosphortransfer_reaction_X = rfr.mol_property_pairs_from_mol_def_and_reaction(mol_def_X,
#                                                                                                  rfr.source_state_set_from_reaction(phosphortransfer_reaction_X))
#
#     print(mol_prop_pair_phosphortransfer_reaction_X)
#
#     assert len(mol_prop_pair_phosphortransfer_reaction_X) == 2
#     assert len(mol_prop_pair_phosphortransfer_reaction_X[0]) == 2
#     assert len(mol_prop_pair_phosphortransfer_reaction_X[1]) == 2
#
#     # !TEST mol_def_E!

@pytest.fixture
def simple_system():
    phosphorylation_reaction_X = rxncon_from_string.reaction_from_string('A_p+_X')
    phosphorylation_reaction_X_at_resi = rxncon_from_string.reaction_from_string('C_p+_X_[(r)]')
    ubiquitination_reaction_X_at_resi = rxncon_from_string.reaction_from_string('D_ub+_X_[(r)]')
    dephosphorylation_reation_X = rxncon_from_string.reaction_from_string('B_p-_X')
    phosphortransfer_reaction_X = rxncon_from_string.reaction_from_string('E_pt_X')

    binding_reaction_A_X = rxncon_from_string.reaction_from_string('A_ppi_X_[d]')
    binding_state_A_X = rxncon_from_string.state_from_string('A--X')

    binding_reaction_Y_X = rxncon_from_string.reaction_from_string('Y_ppi_X_[d]')
    phosphorylated_state = rxncon_from_string.state_from_string('X-{p}')

    binding_contingency1 = con.Contingency(binding_reaction_Y_X,
                                           con.ContingencyType.requirement,
                                           eff.StateEffector(phosphorylated_state))

    binding_contingency2 = con.Contingency(binding_reaction_Y_X,
                                           con.ContingencyType.inhibition,
                                           eff.StateEffector(binding_state_A_X))  # X_ppi_Y; x A--X

    reactions = [phosphorylation_reaction_X, phosphorylation_reaction_X_at_resi, ubiquitination_reaction_X_at_resi,
                 dephosphorylation_reation_X, phosphortransfer_reaction_X, binding_reaction_A_X, binding_reaction_Y_X]

    contingencies = [binding_contingency1, binding_contingency2]

    return RxnConSystem(reactions, contingencies)


@pytest.fixture
def rxn_systems():

    return [
        RxnConSystem([reaction_from_string('A_p+_X')], []),
        RxnConSystem([reaction_from_string('C_p+_X_[(r)]')], []),
        RxnConSystem([reaction_from_string('D_ub+_X_[(r)]'),
                      reaction_from_string('C_p+_X_[(r)]')], []),
        RxnConSystem([reaction_from_string('B_p-_X')], []),
        RxnConSystem([reaction_from_string('E_pt_X')], []),
        RxnConSystem([reaction_from_string('A_ppi_X_[d]'),
                      reaction_from_string('Y_ppi_X_[d]')], []),
        RxnConSystem([reaction_from_string('A_ppi_X_[d]'),
                      reaction_from_string('A_p+_X')],
                     [con.Contingency(reaction_from_string('A_p+_X'),
                                      con.ContingencyType.requirement,
                                      eff.StateEffector(state_from_string('X-{p}')))
                      ]),
        RxnConSystem([reaction_from_string('Y_ppi_X_[d]'),
                      reaction_from_string('Z_ppi_X'),
                      reaction_from_string('A_p+_X')],
                     [con.Contingency(reaction_from_string('Y_ppi_X_[d]'),
                                      con.ContingencyType.requirement,
                                      eff.StateEffector(state_from_string('X-{p}'))),
                      con.Contingency(reaction_from_string('Y_ppi_X_[d]'),
                                      con.ContingencyType.requirement,
                                      eff.StateEffector(state_from_string('Z--X')))])
    ]


@pytest.fixture
def expected_rules(A_pplus_X_expected_rule_system, C_pplus_X_residue_rule_system, D_ubplus_X_residue_C_pplus_X_residue_rule_system,
                   B_pminus_X_expected_rule_system, E_pt_X_expected_rule_system, A_ppi_X_d_and_Y_ppi_X_d_expected_rule_system,
                   A_ppi_X_d_and_A_pplus_X_and_A_ppi_X_d_contingencies_X_pplus_expected_rule_system,
                   Y_ppi_X_d_and_Z_ppi_X_and_Y_ppi_X_d_contingencies_X_pplus_and_Z_bound_X_expected_rule_system):
    return [A_pplus_X_expected_rule_system,
            C_pplus_X_residue_rule_system,
            D_ubplus_X_residue_C_pplus_X_residue_rule_system,
            B_pminus_X_expected_rule_system,
            E_pt_X_expected_rule_system,
            A_ppi_X_d_and_Y_ppi_X_d_expected_rule_system,
            A_ppi_X_d_and_A_pplus_X_and_A_ppi_X_d_contingencies_X_pplus_expected_rule_system,
            Y_ppi_X_d_and_Z_ppi_X_and_Y_ppi_X_d_contingencies_X_pplus_and_Z_bound_X_expected_rule_system
            ]


@pytest.fixture
def Y_ppi_X_d_and_Z_ppi_X_and_Y_ppi_X_d_contingencies_X_pplus_and_Z_bound_X_expected_rule_system(Y_ppi_X_d_reaction, Z_ppi_X_reaction, A_pplus_X_reaction,
                                                                                                        Y_ppi_X_d_contingencies_X_pplus_and_Z_bound_X):
    molecule_definition = mdr.MoleculeDefinitionSupervisor(
        RxnConSystem([Y_ppi_X_d_reaction, Z_ppi_X_reaction, A_pplus_X_reaction],
                     Y_ppi_X_d_contingencies_X_pplus_and_Z_bound_X)).molecule_definitions

    association_def_Y = list(molecule_definition['Y'].association_defs)[0]
    association_def_Z = list(molecule_definition['Z'].association_defs)[0]
    assoc_def_X_bound_Y = [x for x in molecule_definition['X'].association_defs if x.spec.domain == "d"]
    assoc_def_X_bound_Y = assoc_def_X_bound_Y[0]

    assoc_def_X_bound_Z = [x for x in molecule_definition['X'].association_defs if x.spec.domain == "Zassoc"]
    assoc_def_X_bound_Z = assoc_def_X_bound_Z[0]

    modification_def_X = list(molecule_definition['X'].modification_defs)[0]

    association_property_instance_Z_bound_X = moi.AssociationPropertyInstance(association_def_Z, moi.OccupationStatus.occupied_known_partner, assoc_def_X_bound_Z.spec)
    association_property_instance_X_bound_Z = moi.AssociationPropertyInstance(assoc_def_X_bound_Z, moi.OccupationStatus.occupied_known_partner, association_def_Z.spec)

    association_property_instance_X_bound_Y = moi.AssociationPropertyInstance(assoc_def_X_bound_Y, moi.OccupationStatus.occupied_known_partner, association_def_Y.spec)
    association_property_instance_Y_bound_X = moi.AssociationPropertyInstance(association_def_Y, moi.OccupationStatus.occupied_known_partner, assoc_def_X_bound_Y.spec)
                # RULE: Z_PPI_X
    return [rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                set(),
                                                                {moi.AssociationPropertyInstance(assoc_def_X_bound_Z, moi.OccupationStatus.not_occupied, None)},
                                                                None),
                                           ),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['Z'],
                                                                set(), {moi.AssociationPropertyInstance(association_def_Z, moi.OccupationStatus.not_occupied, None)},
                                                                None))
                      ],
                     [rbm.ComplexReactant([moi.MoleculeInstance(molecule_definition['X'],
                                                               set(), {association_property_instance_X_bound_Z}, None),
                                           moi.MoleculeInstance(molecule_definition['Z'],
                                                                set(), {association_property_instance_Z_bound_X}, None)],
                                          [rbm.Binding((0, association_property_instance_X_bound_Z),
                                                       (1, association_property_instance_Z_bound_X))]
                                          )
                      ],
                     rbm.Arrow.reversible,
                     [
                        rbm.Parameter('kf_{0}'.format(str(Z_ppi_X_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(Z_ppi_X_reaction)), None)
                     ]),
            # RULE: A_P+_X
        rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['A'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X,moi.Modifier.unmodified)},
                                                                set(), None))
                      ],  # left_handside
                     [rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['A'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                 {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.phosphorylated)},
                                                                 set(), None))
                      ],  # right_handside
                     rbm.Arrow.reversible,  # arrow_type  #  should be reversible but is unidirectional hence ->
                     [
                        rbm.Parameter('kf_{0}'.format(str(A_pplus_X_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(A_pplus_X_reaction)), None)
                     ]
                     ),

        # RULE: Y_PPI_X_d ! X-{P} A--X
        rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['Y'],
                                                                set(), {moi.AssociationPropertyInstance(association_def_Y, moi.OccupationStatus.not_occupied, None)},
                                                                None)),
                  rbm.ComplexReactant([moi.MoleculeInstance(molecule_definition['X'],
                                                            {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.phosphorylated)},
                                                            {association_property_instance_X_bound_Z,
                                                             moi.AssociationPropertyInstance(assoc_def_X_bound_Y, moi.OccupationStatus.not_occupied, None)}, None),
                                           moi.MoleculeInstance(molecule_definition['Z'],
                                                                set(), {association_property_instance_Z_bound_X}, None)],
                                          [rbm.Binding((0, association_property_instance_X_bound_Z),
                                                       (1, association_property_instance_Z_bound_X))]
                                          )
                      ],
                     [rbm.ComplexReactant([moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.phosphorylated)},
                                                                {association_property_instance_X_bound_Y,
                                                                 association_property_instance_X_bound_Z},
                                                                None),
                                           moi.MoleculeInstance(molecule_definition['Y'],
                                                                set(), {association_property_instance_Y_bound_X},
                                                                None),
                                           moi.MoleculeInstance(molecule_definition['Z'],
                                                                set(), {association_property_instance_Z_bound_X}, None)],
                                          [rbm.Binding((0, association_property_instance_X_bound_Z),  # Zassoc
                                                       (2, association_property_instance_Z_bound_X)), # Xassoc
                                           rbm.Binding((0, association_property_instance_X_bound_Y),  # d
                                                       (1, association_property_instance_Y_bound_X)), # XAssoc
                                           ]),
                      ],
                     rbm.Arrow.reversible,
                     [
                        rbm.Parameter('kf_{0}'.format(str(Y_ppi_X_d_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(Y_ppi_X_d_reaction)), None)
                     ])
        ]


@pytest.fixture
def A_ppi_X_d_and_A_pplus_X_and_A_ppi_X_d_contingencies_X_pplus_expected_rule_system(A_ppi_X_d_reaction, A_pplus_X_reaction, A_ppi_X_d_contingencies_X_pplus):
    molecule_definition = mdr.MoleculeDefinitionSupervisor(
        RxnConSystem([A_ppi_X_d_reaction, A_pplus_X_reaction], A_ppi_X_d_contingencies_X_pplus)).molecule_definitions
    association_def_A = list(molecule_definition['A'].association_defs)[0]
    association_def_X = list(molecule_definition['X'].association_defs)[0]
    modification_def_X = list(molecule_definition['X'].modification_defs)[0]

    association_property_instance_A_bound_X = moi.AssociationPropertyInstance(association_def_A, moi.OccupationStatus.occupied_known_partner, association_def_X.spec)
    association_property_instance_X_bound_A = moi.AssociationPropertyInstance(association_def_X, moi.OccupationStatus.occupied_known_partner, association_def_A.spec)
            # RULE: A_PPI_X_d ! X-{P}
    return [rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['A'],
                                                                set(), {moi.AssociationPropertyInstance(association_def_A, moi.OccupationStatus.not_occupied, None)},
                                                                None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X,moi.Modifier.phosphorylated)},
                                                                {moi.AssociationPropertyInstance(association_def_X, moi.OccupationStatus.not_occupied, None)},
                                                                None))
                      ],
                     [rbm.ComplexReactant([moi.MoleculeInstance(molecule_definition['A'],
                                                                set(), {association_property_instance_A_bound_X}, None),
                                           moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X,moi.Modifier.phosphorylated)},
                                                                {association_property_instance_X_bound_A}, None)],
                                          [rbm.Binding((0, association_property_instance_A_bound_X),
                                                       (1, association_property_instance_X_bound_A))]
                                          ),
                      ],
                     rbm.Arrow.reversible,
                     [
                        rbm.Parameter('kf_{0}'.format(str(A_ppi_X_d_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(A_ppi_X_d_reaction)), None)
                     ]),
            #RULE: A_P+_X
            rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['A'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X,moi.Modifier.unmodified)},
                                                                set(), None))],  # left_handside
                     [rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['A'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                 {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.phosphorylated)},
                                                                 set(), None))],  # right_handside
                     rbm.Arrow.reversible,  # arrow_type  #  should be reversible but is unidirectional hence ->
                     [
                        rbm.Parameter('kf_{0}'.format(str(A_pplus_X_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(A_pplus_X_reaction)), None)
                     ]
                     )
            ]


@pytest.fixture
def A_ppi_X_d_and_Y_ppi_X_d_expected_rule_system(A_ppi_X_d_reaction, Y_ppi_X_d_reaction):
    molecule_definition = mdr.MoleculeDefinitionSupervisor(RxnConSystem([A_ppi_X_d_reaction, Y_ppi_X_d_reaction], [])).molecule_definitions
    association_def_A = list(molecule_definition['A'].association_defs)[0]
    association_def_Y = list(molecule_definition['Y'].association_defs)[0]
    association_def_X = list(molecule_definition['X'].association_defs)[0]

    association_property_instance_A_bound_X = moi.AssociationPropertyInstance(association_def_A, moi.OccupationStatus.occupied_known_partner, association_def_X.spec)
    association_property_instance_X_bound_A = moi.AssociationPropertyInstance(association_def_X, moi.OccupationStatus.occupied_known_partner, association_def_A.spec)
    #
    association_property_instance_X_bound_Y = moi.AssociationPropertyInstance(association_def_X, moi.OccupationStatus.occupied_known_partner, association_def_Y.spec)
    association_property_instance_Y_bound_X = moi.AssociationPropertyInstance(association_def_Y, moi.OccupationStatus.occupied_known_partner, association_def_X.spec)
            # RULE: A_PPI_X_d
    return [rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['A'],
                                                                set(), {moi.AssociationPropertyInstance(association_def_A, moi.OccupationStatus.not_occupied, None)},
                                                                None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                set(), {moi.AssociationPropertyInstance(association_def_X, moi.OccupationStatus.not_occupied, None)},
                                                                None))
                      ],
                     [rbm.ComplexReactant([moi.MoleculeInstance(molecule_definition['A'],
                                                                set(), {association_property_instance_A_bound_X}, None),
                                           moi.MoleculeInstance(molecule_definition['X'],
                                                                set(), {association_property_instance_X_bound_A}, None)],
                                          [rbm.Binding((0, association_property_instance_A_bound_X),
                                                       (1, association_property_instance_X_bound_A))]
                                          ),
                      ],
                     rbm.Arrow.reversible,
                     [
                        rbm.Parameter('kf_{0}'.format(str(A_ppi_X_d_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(A_ppi_X_d_reaction)), None)
                     ]),
            # RULE: Y_PPI_X_d
            rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                set(), {moi.AssociationPropertyInstance(association_def_X, moi.OccupationStatus.not_occupied, None)},
                                                                None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['Y'],
                                                                set(), {moi.AssociationPropertyInstance(association_def_Y, moi.OccupationStatus.not_occupied, None)},
                                                                None))
                      ],
                     [rbm.ComplexReactant([moi.MoleculeInstance(molecule_definition['X'],
                                                                set(), {association_property_instance_X_bound_Y},
                                                                None),
                                           moi.MoleculeInstance(molecule_definition['Y'],
                                                                set(), {association_property_instance_Y_bound_X},
                                                                None)],
                                          [rbm.Binding((0, association_property_instance_X_bound_Y),
                                                       (1, association_property_instance_Y_bound_X))]),
                      ],
                     rbm.Arrow.reversible,
                     [
                        rbm.Parameter('kf_{0}'.format(str(Y_ppi_X_d_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(Y_ppi_X_d_reaction)), None)
                     ])
            ]


@pytest.fixture
def E_pt_X_expected_rule_system(E_pt_X_reaction):
    # todo: B_pt_E are two reactions in one E_p+_X -> X_[Eside] and X_p-_E -> E_[Xside]
    # todo: B_[n]_apt_B_[m] auto phosphortransfer B is the same molecule B_[n]_p+_B_[m] -> B_[m] and B_[m]_p-_B_[n] -> B_B[n]
    # todo: B_apt_B auto phosphortransfer B is the same molecule B_p+_B -> B_[Bsite1] and B_p-_B -> B_B[Site2]

    molecule_definition = mdr.MoleculeDefinitionSupervisor(RxnConSystem([E_pt_X_reaction], [])).molecule_definitions
    modification_def_E = list(molecule_definition['E'].modification_defs)[0]
    modification_def_X = list(molecule_definition['X'].modification_defs)[0]
            # RULE: E_PT_X
    return [rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['E'],
                                                                {moi.ModificationPropertyInstance(modification_def_E, moi.Modifier.phosphorylated)}
                                                                , set(), None)),
                     rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                               {moi.ModificationPropertyInstance(modification_def_X,moi.Modifier.unmodified)},
                                                               set(), None))],  # left_reactant
                     [rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['E'],
                                                                {moi.ModificationPropertyInstance(modification_def_E, moi.Modifier.unmodified)},
                                                                set(), None)),
                     rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                               {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.phosphorylated)},
                                                               set(), None))],  # right_reactant
                     rbm.Arrow.reversible,  # arrow_type  #  should be reversible but is unidirectional hence ->
                     rfr.parameters_from_reaction_and_quant_conts(E_pt_X_reaction, [])
                      )
            ]


@pytest.fixture
def B_pminus_X_expected_rule_system(B_pminus_X_reaction):

    molecule_definition = mdr.MoleculeDefinitionSupervisor(RxnConSystem([B_pminus_X_reaction], [])).molecule_definitions
    modification_def_X = list(molecule_definition['X'].modification_defs)[0]
            # RULE: B_P-_X
    return [rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['B'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X,moi.Modifier.phosphorylated)},
                                                                set(), None))],  # left_reactant
                     [rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['B'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                 {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.unmodified)},
                                                                 set(), None))],  # right_reactant
                     rbm.Arrow.reversible,  # arrow_type  #  should be reversible but is unidirectional hence ->

                     [
                        rbm.Parameter('kf_{0}'.format(str(B_pminus_X_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(B_pminus_X_reaction)), None)
                     ]
                     )
            ]


@pytest.fixture
def A_pplus_X_expected_rule_system(A_pplus_X_reaction):

    molecule_definition = mdr.MoleculeDefinitionSupervisor(RxnConSystem([A_pplus_X_reaction], [])).molecule_definitions
    modification_def_X = list(molecule_definition['X'].modification_defs)[0]

            # RULE: A_P+_X
    return [rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['A'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X,moi.Modifier.unmodified)},
                                                                set(), None))],  # left_reactant
                     [rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['A'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                 {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.phosphorylated)},
                                                                 set(), None))],  # right_reactant
                     rbm.Arrow.reversible,  # arrow_type  #  should be reversible but is unidirectional hence ->
                     [
                        rbm.Parameter('kf_{0}'.format(str(A_pplus_X_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(A_pplus_X_reaction)), None)
                     ]
                     )
            ]


@pytest.fixture
def C_pplus_X_residue_rule_system(C_pplus_X_residue_reaction):

    molecule_definition = mdr.MoleculeDefinitionSupervisor(RxnConSystem([C_pplus_X_residue_reaction], [])).molecule_definitions
    modification_def_X = list(molecule_definition['X'].modification_defs)[0]
            # RULE: C_P+_X
    return [rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['C'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.unmodified)},
                                                                set(), None))],  # left_reactant
                     [rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['C'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.phosphorylated)},
                                                                set(), None))],  # right_reactant
                     rbm.Arrow.reversible,  # arrow_type  #  should be reversible but is unidirectional hence ->
                     [
                        rbm.Parameter('kf_{0}'.format(str(C_pplus_X_residue_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(C_pplus_X_residue_reaction)), None)
                     ]
                     )
            ]


@pytest.fixture
def D_ubplus_X_residue_C_pplus_X_residue_rule_system(D_ubplus_X_residue_reaction, C_pplus_X_residue_reaction):

    molecule_definition = mdr.MoleculeDefinitionSupervisor(
        RxnConSystem([D_ubplus_X_residue_reaction, C_pplus_X_residue_reaction], [])).molecule_definitions
    assert len(molecule_definition['X'].modification_defs) == 1  # there should be only one residue which is modified by two reactions
    modification_def_X = list(molecule_definition['X'].modification_defs)[0]
            #RULE: D_UB+_X_(r)
    return [rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['C'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.unmodified)},
                                                                set(), None))],  # left_reactant
                     [rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['C'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.phosphorylated)},
                                                                set(), None))],  # right_reactant
                     rbm.Arrow.reversible,  # arrow_type  #  should be reversible but is unidirectional hence ->
                     [
                        rbm.Parameter('kf_{0}'.format(str(C_pplus_X_residue_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(C_pplus_X_residue_reaction)), None)
                     ]
                      ),
            # RULE: C_P+_X_(r)
            rbm.Rule([rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['D'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.unmodified)},
                                                                set(), None))],  # left_reactant
                     [rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['D'], set(), set(), None)),
                      rbm.MoleculeReactant(moi.MoleculeInstance(molecule_definition['X'],
                                                                {moi.ModificationPropertyInstance(modification_def_X, moi.Modifier.ubiquitinated)},
                                                                set(), None))],  # right_reactant
                     rbm.Arrow.reversible,  # arrow_type  #  should be reversible but is unidirectional hence ->
                     [
                        rbm.Parameter('kf_{0}'.format(str(D_ubplus_X_residue_reaction)), None),
                        rbm.Parameter('kr_{0}'.format(str(D_ubplus_X_residue_reaction)), None)
                     ]
                     )
            ]
