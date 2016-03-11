import pytest
import typing as tp
import rxncon.core.rxncon_system as rxs
import rxncon.simulation
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.core.contingency as con
import rxncon.core.specification as spec
import rxncon.core.effector as eff
import rxncon.simulation.rule_based.rbm_from_rxncon as rfr
from rxncon.core import contingency as con, effector, effector, effector, effector, rxncon_system
from rxncon.input.quick import quick as qui
from rxncon.syntax import rxncon_from_string as rfs
from rxncon.venntastic import sets as venn
import rxncon.semantics.molecule_definition_from_rxncon as mdr
import rxncon.semantics.molecule_instance as moi
import rxncon.semantics.molecule_instance_from_rxncon as mfr

# MASTERTEST testing the lhs/rhs MoleculeInstances that appear.

def test_mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(simple_system: rxs.RxnConSystem):

    phosphorylation_reaction_X = simple_system.reactions[0]

    mol_defs = mdr.MoleculeDefinitionSupervisor(simple_system).molecule_definitions
    mol_def_X = mol_defs['X']

    strict_conts = simple_system.strict_contingencies_for_reaction(phosphorylation_reaction_X)

    phosphorylation_reaction_X_pairs = rfr.mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(mol_def_X, phosphorylation_reaction_X, strict_conts)
    mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "Asite"]

    expected_lhs = moi.MoleculeInstance(mol_def_X, {moi.ModificationPropertyInstance(mod_defs[0], moi.Modifier("u"))}, set(), None)
    expected_rhs = moi.MoleculeInstance(mol_def_X, {moi.ModificationPropertyInstance(mod_defs[0], moi.Modifier("p"))}, set(), None)
    assert phosphorylation_reaction_X_pairs == [(expected_lhs, expected_rhs)]

    phosphorylation_reaction_X_at_resi = simple_system.reactions[1]
    mol_defs = mdr.MoleculeDefinitionSupervisor(simple_system).molecule_definitions
    mol_def_X = mol_defs['X']
    strict_conts = simple_system.strict_contingencies_for_reaction(phosphorylation_reaction_X_at_resi)

    phosphorylation_reaction_X_at_resi_pairs = rfr.mol_instance_pairs_from_mol_def_and_reaction_and_contingencies(mol_def_X, phosphorylation_reaction_X_at_resi, strict_conts)

    mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "r"]

    expected_lhs_phosphorylation_reaction_X_at_resi = moi.MoleculeInstance(mol_def_X, {moi.ModificationPropertyInstance(mod_defs[0], moi.Modifier("u"))}, set(), None)
    expected_rhs_phosphorylation_reaction_X_at_resi = moi.MoleculeInstance(mol_def_X, {moi.ModificationPropertyInstance(mod_defs[0], moi.Modifier("p"))}, set(), None)
    assert phosphorylation_reaction_X_at_resi_pairs == [(expected_lhs_phosphorylation_reaction_X_at_resi,
                                                         expected_rhs_phosphorylation_reaction_X_at_resi)]

def test_set_of_states_from_single_state_effector(simple_system: rxs.RxnConSystem):
    cont1 = simple_system.contingencies[0]
    cont2 = simple_system.contingencies[1]

    state_set_of_effector_cont1 = rfr.state_set_from_effector(cont1.effector)
    state_set_of_effector_cont2 = rfr.state_set_from_effector(cont2.effector)

    assert cont1.target == rfs.reaction_from_string('Y_ppi_X_[d]')
    assert cont1.type == con.ContingencyType.requirement
    assert state_set_of_effector_cont1.is_equivalent_to(venn.PropertySet(rfs.state_from_string('X-{P}')))

    assert cont2.target == rfs.reaction_from_string('Y_ppi_X_[d]')
    assert cont2.type == con.ContingencyType.inhibition
    assert state_set_of_effector_cont2.is_equivalent_to(venn.PropertySet(rfs.state_from_string('A--X')))


def test_set_of_states_from_nested_AND_effectors():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)

    assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.PropertySet(expected_b_dash_f)))


def test_set_of_states_from_nested_OR_effectors():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR A--C
                        <comp>; OR C--E
                        <comp>; OR B--F""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Union(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.PropertySet(expected_b_dash_f)))


def test_set_of_states_from_nested_AND_OR_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND <c1>
                        <comp>; AND <c2>
                        <c1>; OR A--C
                        <c1>; OR C--E
                        <c2>; OR B--F
                        <c2>; OR B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Intersection(venn.Union(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Union(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_nested_OR_AND_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; AND B--F
                        <c2>; AND B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")


    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Intersection(venn.PropertySet(expected_b_dash_f),venn.PropertySet(expected_b_dash_d))))


def test_set_of_states_from_complex_effector():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; OR <c1>
                        <comp>; OR <c2>
                        <c1>; AND A--C
                        <c1>; AND C--E
                        <c2>; NOT <c3>
                        <c3>; AND B--F
                        <c3>; AND B--D""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")
    expected_b_dash_d = rfs.state_from_string("B--D")

    rxncon = quick.rxncon_system

    set_of_state_AND_effector = rfr.state_set_from_effector(rxncon.contingencies[0].effector)
    assert set_of_state_AND_effector.is_equivalent_to(venn.Union(venn.Intersection(venn.PropertySet(expected_a_dash_c),
                                                                                          venn.PropertySet(expected_c_dash_e)),
                                                                        venn.Complement(venn.Intersection(venn.PropertySet(expected_b_dash_d),
                                                                                                          venn.PropertySet(expected_b_dash_f)))))


def test_state_set_from_contingencies(simple_system):


    a_ppi_b = rfs.reaction_from_string('A_ppi_B')
    a_dash_b = rfs.state_from_string('A--B')

    a_ppi_c = rfs.reaction_from_string('A_ppi_C')
    a_dash_c = rfs.state_from_string('A--C')

    b_ppi_e = rfs.reaction_from_string('B_ppi_E')

    b_pplus_e = rfs.reaction_from_string('B_p+_E')
    e_pplus = rfs.state_from_string("E-{P}")

    cont_b_dash_e = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(a_dash_b))  # B_ppi_E; ! A--B
    cont_e_pplus = con.Contingency(b_ppi_e, con.ContingencyType.requirement, eff.StateEffector(e_pplus))  # B_ppi_E; ! E-{P}
    cont_a_ppi_b = con.Contingency(a_ppi_b, con.ContingencyType.inhibition, eff.StateEffector(a_dash_c))  # A_ppi_B; x A--C

    rxncon = rxs.RxnConSystem([b_ppi_e, a_ppi_b, a_ppi_c, b_pplus_e], [cont_e_pplus, cont_b_dash_e, cont_a_ppi_b])

    strict_contingencies_state_set_b_ppi_e = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(b_ppi_e))
    strict_contingencies_state_set_a_ppi_b = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_b))
    strict_contingencies_state_set_a_ppi_c = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(a_ppi_c))
    strict_contingencies_state_set_b_pplus_e = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(b_pplus_e))

    expected_b_ppi_e_strict_cont = venn.Intersection(venn.PropertySet(e_pplus), venn.PropertySet(a_dash_b))
    assert strict_contingencies_state_set_b_ppi_e.is_equivalent_to(expected_b_ppi_e_strict_cont)

    expected_a_ppi_b_strict_cont = venn.Complement(venn.PropertySet(a_dash_c))
    assert strict_contingencies_state_set_a_ppi_b.is_equivalent_to(expected_a_ppi_b_strict_cont)

    expected_a_ppi_c_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_a_ppi_c.is_equivalent_to(expected_a_ppi_c_strict_cont)

    expected_b_pplus_e_strict_cont = venn.UniversalSet()
    assert strict_contingencies_state_set_b_pplus_e.is_equivalent_to(expected_b_pplus_e_strict_cont)


def test_set_of_states_from_contingencies_FOR_quant():
    # todo
    pass


def test_set_of_states_from_strict_contingencies_AND():
    quick = qui.Quick("""A_ppi_B; ! <comp>
                        <comp>; AND A--C
                        <comp>; AND C--E
                        <comp>; AND B--F""")

    expected_a_dash_c = rfs.state_from_string("A--C")
    expected_c_dash_e = rfs.state_from_string("C--E")
    expected_b_dash_f = rfs.state_from_string("B--F")

    rxncon = quick.rxncon_system

    strict_cont_state_set = rfr.state_set_from_contingencies(rxncon.strict_contingencies_for_reaction(rfs.reaction_from_string('A_ppi_B')))

    assert strict_cont_state_set.is_equivalent_to(venn.Intersection(venn.Intersection(venn.PropertySet(expected_a_dash_c), venn.PropertySet(expected_c_dash_e)),
                                                                    venn.PropertySet(expected_b_dash_f)))


def test_source_state_set_from_reaction(simple_system: rxs.RxnConSystem):

    phosphorylation_reaction_X = simple_system.reactions[0]
    phosphorylation_reaction_X_at_resi = simple_system.reactions[1]
    ubiquitination_reaction_X_at_resi = simple_system.reactions[2]
    dephosphorylation_reaction_X = simple_system.reactions[3]
    phosphortransfer_reaction_X = simple_system.reactions[4]
    binding_reaction_A_X = simple_system.reactions[5]
    binding_reaction_Y_X = simple_system.reactions[6]

    assert rfr.source_state_set_from_reaction(phosphorylation_reaction_X).is_equivalent_to(venn.Complement(
                                                                                                venn.PropertySet(
                                                                                                    rfs.state_from_string(
                                                                                                        'X-{P}'))))
    assert rfr.source_state_set_from_reaction(phosphorylation_reaction_X_at_resi).is_equivalent_to(venn.Complement(
                                                                                                        venn.PropertySet(
                                                                                                            rfs.state_from_string('X_(r)-{P}'))))
    assert rfr.source_state_set_from_reaction(ubiquitination_reaction_X_at_resi).is_equivalent_to(venn.Complement(
                                                                                                        venn.PropertySet(
                                                                                                            rfs.state_from_string('X_(r)-{Ub}'))))
    assert rfr.source_state_set_from_reaction(dephosphorylation_reaction_X).is_equivalent_to(venn.PropertySet(rfs.state_from_string('X-{P}')))

    # todo: B_pt_E are two reactions in one E_p+_X -> X_[Eside] and X_p-_E -> E_[Xside]
    # todo: B_[n]_apt_B_[m] auto phosphortransfer B is the same molecule B_[n]_p+_B_[m] -> B_[m] and B_[m]_p-_B_[n] -> B_B[n]
    # todo: B_apt_B auto phosphortransfer B is the same molecule B_p+_B -> B_[Bsite1] and B_p-_B -> B_B[Site2]

    assert rfr.source_state_set_from_reaction(phosphortransfer_reaction_X).is_equivalent_to(venn.Intersection(venn.Complement(venn.PropertySet(rfs.state_from_string('X-{p}'))),
                                                                                                              venn.PropertySet(rfs.state_from_string('E-{p}'))))
    assert rfr.source_state_set_from_reaction(binding_reaction_A_X).is_equivalent_to(venn.Complement(venn.PropertySet(rfs.state_from_string('A--X_[d]'))))
    assert rfr.source_state_set_from_reaction(binding_reaction_Y_X).is_equivalent_to(venn.Complement(venn.PropertySet(rfs.state_from_string('Y--X_[d]'))))


def test_mol_property_pairs_from_mol_def_and_source_state_set(simple_system):
    mol_defs = mdr.MoleculeDefinitionSupervisor(simple_system).molecule_definitions
    mol_def_X = mol_defs['X']

    phosphorylation_reaction_X = simple_system.reactions[0]
    phosphorylation_reaction_X_at_resi = simple_system.reactions[1]
    ubiquitination_reaction_X_at_resi = simple_system.reactions[2]
    dephosphorylation_reaction_X = simple_system.reactions[3]
    phosphortransfer_reaction_X = simple_system.reactions[4]
    binding_reaction_A_X = simple_system.reactions[5]
    binding_reaction_Y_X = simple_system.reactions[6]


    mol_prop_pair_phosphorylation_reaction_X = rfr.mol_property_pairs_from_mol_def_and_source_state_set(mol_def_X,
                                                                                                        rfr.source_state_set_from_reaction(phosphorylation_reaction_X))
    assert len(mol_prop_pair_phosphorylation_reaction_X) == 1
    assert len(mol_prop_pair_phosphorylation_reaction_X[0]) == 2

    mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "Asite"]
    assert len(mod_defs) == 1
    mod_def_Asite = mod_defs[0]
    assert mol_prop_pair_phosphorylation_reaction_X[0][0] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_Asite, moi.Modifier.unmodified))
    assert mol_prop_pair_phosphorylation_reaction_X[0][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_Asite, moi.Modifier.phosphorylated))


    mol_prop_pair_phosphorylation_reaction_X_at_resi = rfr.mol_property_pairs_from_mol_def_and_source_state_set(mol_def_X,
                                                                                                                rfr.source_state_set_from_reaction(phosphorylation_reaction_X_at_resi))

    assert len(mol_prop_pair_phosphorylation_reaction_X_at_resi) == 2
    assert len(mol_prop_pair_phosphorylation_reaction_X_at_resi[0]) == 2
    assert len(mol_prop_pair_phosphorylation_reaction_X_at_resi[1]) == 2

    mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "r"]
    assert len(mod_defs) == 1
    mod_def_r = mod_defs[0]

    assert mol_prop_pair_phosphorylation_reaction_X_at_resi[0][0] != mol_prop_pair_phosphorylation_reaction_X_at_resi[1][0]
    assert mol_prop_pair_phosphorylation_reaction_X_at_resi[0][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_r, moi.Modifier.phosphorylated))
    assert mol_prop_pair_phosphorylation_reaction_X_at_resi[1][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_r, moi.Modifier.phosphorylated))


    mol_prop_pair_ubiquitination_reaction_X_at_resi = rfr.mol_property_pairs_from_mol_def_and_source_state_set(mol_def_X,
                                                                                                               rfr.source_state_set_from_reaction(ubiquitination_reaction_X_at_resi))

    assert len(mol_prop_pair_ubiquitination_reaction_X_at_resi) == 2
    assert len(mol_prop_pair_ubiquitination_reaction_X_at_resi[0]) == 2
    assert len(mol_prop_pair_ubiquitination_reaction_X_at_resi[1]) == 2


    assert mol_prop_pair_ubiquitination_reaction_X_at_resi[0][0] != mol_prop_pair_ubiquitination_reaction_X_at_resi[1][0]
    assert mol_prop_pair_ubiquitination_reaction_X_at_resi[0][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_r, moi.Modifier.ubiquitinated))
    assert mol_prop_pair_ubiquitination_reaction_X_at_resi[1][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_r, moi.Modifier.ubiquitinated))


    mol_prop_pair_dephosphorylation_reaction_X = rfr.mol_property_pairs_from_mol_def_and_source_state_set(mol_def_X,
                                                                                                          rfr.source_state_set_from_reaction(dephosphorylation_reaction_X))


    assert len(mol_prop_pair_dephosphorylation_reaction_X) == 1
    assert len(mol_prop_pair_dephosphorylation_reaction_X[0]) == 2

    mod_defs = [x for x in mol_def_X.modification_defs if x.spec.residue == "Bsite"]
    assert len(mod_defs) == 1
    mod_def_Bsite = mod_defs[0]

    assert mol_prop_pair_dephosphorylation_reaction_X[0][0] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_Bsite, moi.Modifier.phosphorylated))
    assert mol_prop_pair_dephosphorylation_reaction_X[0][1] == venn.PropertySet(moi.ModificationPropertyInstance(mod_def_Bsite, moi.Modifier.unmodified))


    mol_prop_pair_binding_reaction_A_X = rfr.mol_property_pairs_from_mol_def_and_source_state_set(mol_def_X,
                                                                                                  rfr.source_state_set_from_reaction(binding_reaction_A_X))

    assert len(mol_prop_pair_binding_reaction_A_X) == 2
    assert len(mol_prop_pair_binding_reaction_A_X[0]) == 2
    assert len(mol_prop_pair_binding_reaction_A_X[1]) == 2

    assoc_defs = [x for x in mol_def_X.association_defs if x.spec.domain == "d"]
    assert len(assoc_defs) == 1
    assoc_def_d = assoc_defs[0]

    assert mol_prop_pair_binding_reaction_A_X[0][0] != mol_prop_pair_binding_reaction_A_X[1][0]
    assert mol_prop_pair_binding_reaction_A_X[0][1] == venn.PropertySet(moi.AssociationPropertyInstance(assoc_def_d,
                                                                                                        moi.OccupationStatus.occupied_known_partner,
                                                                                                        spec.Specification('A', 'Xassoc', None, None)))
    assert mol_prop_pair_binding_reaction_A_X[1][1] == venn.PropertySet(moi.AssociationPropertyInstance(assoc_def_d,
                                                                                                        moi.OccupationStatus.occupied_known_partner,
                                                                                                        spec.Specification('A', 'Xassoc', None, None)))
    mol_prop_pair_binding_reaction_Y_X = rfr.mol_property_pairs_from_mol_def_and_source_state_set(mol_def_X,
                                                                                                  rfr.source_state_set_from_reaction(binding_reaction_Y_X))

    assert len(mol_prop_pair_binding_reaction_Y_X) == 2
    assert len(mol_prop_pair_binding_reaction_Y_X[0]) == 2
    assert len(mol_prop_pair_binding_reaction_Y_X[1]) == 2

    assert mol_prop_pair_binding_reaction_Y_X[0][0] != mol_prop_pair_binding_reaction_Y_X[1][0]
    assert mol_prop_pair_binding_reaction_Y_X[0][1] == venn.PropertySet(moi.AssociationPropertyInstance(assoc_def_d,
                                                                                                        moi.OccupationStatus.occupied_known_partner,
                                                                                                        spec.Specification('Y', 'Xassoc', None, None)))
    assert mol_prop_pair_binding_reaction_Y_X[1][1] == venn.PropertySet(moi.AssociationPropertyInstance(assoc_def_d,
                                                                                                        moi.OccupationStatus.occupied_known_partner,
                                                                                                        spec.Specification('Y', 'Xassoc', None, None)))
    #is_property_pair_valid_for_reaction(mol_def, x, reaction)

    mol_prop_pair_phosphortransfer_reaction_X = rfr.mol_property_pairs_from_mol_def_and_source_state_set(mol_def_X,
                                                                                                         rfr.source_state_set_from_reaction(phosphortransfer_reaction_X))

    print(mol_prop_pair_phosphortransfer_reaction_X)

    assert len(mol_prop_pair_phosphortransfer_reaction_X) == 2
    assert len(mol_prop_pair_phosphortransfer_reaction_X[0]) == 2
    assert len(mol_prop_pair_phosphortransfer_reaction_X[1]) == 2

    # !TEST mol_def_E!

@pytest.fixture
def simple_system():
    phosphorylation_reaction_X = rfs.reaction_from_string('A_p+_X')
    phosphorylation_reaction_X_at_resi = rfs.reaction_from_string('C_p+_X_[(r)]')
    ubiquitination_reaction_X_at_resi = rfs.reaction_from_string('D_ub+_X_[(r)]')
    dephosphorylation_reation_X = rfs.reaction_from_string('B_p-_X')
    phosphortransfer_reaction_X = rfs.reaction_from_string('E_pt_X')

    binding_reaction_A_X = rfs.reaction_from_string('A_ppi_X_[d]')
    binding_state_A_X = rfs.state_from_string('A--X')

    binding_reaction_Y_X = rfs.reaction_from_string('Y_ppi_X_[d]')
    phosphorylated_state = rfs.state_from_string('X-{p}')

    binding_contingency1 = con.Contingency(binding_reaction_Y_X,
                                           con.ContingencyType.requirement,
                                           eff.StateEffector(phosphorylated_state))

    binding_contingency2 = con.Contingency(binding_reaction_Y_X,
                                           con.ContingencyType.inhibition,
                                           eff.StateEffector(binding_state_A_X))  # X_ppi_Y; x A--X

    reactions = [phosphorylation_reaction_X, phosphorylation_reaction_X_at_resi, ubiquitination_reaction_X_at_resi,
                 dephosphorylation_reation_X, phosphortransfer_reaction_X, binding_reaction_A_X, binding_reaction_Y_X]

    contingencies = [binding_contingency1, binding_contingency2]

    return rxs.RxnConSystem(reactions, contingencies)
