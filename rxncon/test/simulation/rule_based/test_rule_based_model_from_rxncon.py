
import rxncon.core.rxncon_system as rxs
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.simulation.rule_based.rule_based_model as rbm
import rxncon.simulation.rule_based.rule_based_model_from_rxncon as rfr
import rxncon.semantics.state_flow as flo
import rxncon.input.quick.quick as qui
import rxncon.simulation.rule_based.bngl_export as bng

### MOLECULE DEFINITIONS FROM SINGLE REACTION SYSTEMS ###
def test_single_ppi_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_ppi_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    #mol_def_a = [x for x in mol_defs if x.name == 'A'][0]
    #mol_def_b = [x for x in mol_defs if x.name == 'B'][0]

    assert not mol_defs["A"].modification_defs
    assert not mol_defs["A"].localization_def

    assert len(mol_defs["A"].association_defs) == 1
    assert mol_defs["A"].association_defs[0].domain == 'AssocB'

    assert not mol_defs["B"].modification_defs
    assert not mol_defs["B"].localization_def

    assert len(mol_defs["B"].association_defs) == 1
    assert mol_defs["B"].association_defs[0].domain == 'AssocA'

def test_double_ppi_reaction_mol_defs():
    reaction1 = rfs.reaction_from_string('A_ppi_B')
    reaction2 = rfs.reaction_from_string('A_ppi_C')
    rxncon_sys = rxs.RxnConSystem([reaction1, reaction2], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)


    assert not mol_defs["A"].modification_defs
    assert not mol_defs["A"].localization_def

    assert len(mol_defs["A"].association_defs) == 2
    assert mol_defs["A"].association_defs[0].domain != mol_defs["A"].association_defs[1].domain
    assert mol_defs["A"].association_defs[0].domain in ['AssocB', 'AssocC']
    assert mol_defs["A"].association_defs[1].domain in ['AssocB', 'AssocC']

    assert not mol_defs["B"].modification_defs
    assert not mol_defs["B"].localization_def

    assert len(mol_defs["B"].association_defs) == 1
    assert mol_defs["B"].association_defs[0].domain == 'AssocA'

    assert not mol_defs["C"].modification_defs
    assert not mol_defs["C"].localization_def

    assert len(mol_defs["C"].association_defs) == 1
    assert mol_defs["C"].association_defs[0].domain == 'AssocA'


def test_single_ipi_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_[n]_ipi_A_[d]')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    assert not mol_defs["A"].modification_defs
    assert not mol_defs["A"].localization_def

    assert len(mol_defs["A"].association_defs) == 2
    assert mol_defs["A"].association_defs[0] != mol_defs["A"].association_defs[1]
    assert mol_defs["A"].association_defs[0].domain in ["n", "d"]
    assert mol_defs["A"].association_defs[1].domain in ["n", "d"]


def test_single_modification_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_p+_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)


    assert not mol_defs["A"].association_defs
    assert not mol_defs["A"].modification_defs
    assert not mol_defs["A"].localization_def

    assert not mol_defs["B"].association_defs
    assert not mol_defs["B"].localization_def

    assert len(mol_defs["B"].modification_defs) == 1
    assert mol_defs["B"].modification_defs[0].domain == "ModA"
    assert mol_defs["B"].modification_defs[0].valid_modifiers == ["u","p"]

def test_double_modification_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_p+_B')
    reaction1 = rfs.reaction_from_string('C_p+_B_[c]')
    rxncon_sys = rxs.RxnConSystem([reaction, reaction1], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)


    assert not mol_defs["A"].association_defs
    assert not mol_defs["A"].modification_defs
    assert not mol_defs["A"].localization_def

    assert not mol_defs["C"].association_defs
    assert not mol_defs["C"].modification_defs
    assert not mol_defs["C"].localization_def

    assert not mol_defs["B"].association_defs
    assert not mol_defs["B"].localization_def

    assert len(mol_defs["B"].modification_defs) == 2
    assert mol_defs["B"].modification_defs[0].domain != mol_defs["B"].modification_defs[1].domain
    assert mol_defs["B"].modification_defs[0].domain in ["ModA","c"]
    assert mol_defs["B"].modification_defs[1].domain in ["ModA","c"]

    assert mol_defs["B"].modification_defs[0].valid_modifiers == ["u","p"]
    assert mol_defs["B"].modification_defs[1].valid_modifiers == ["u","p"]



def test_single_synth_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_syn_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    assert not mol_defs["A"].modification_defs
    assert not mol_defs["A"].localization_def
    assert not mol_defs["A"].association_defs

    assert not mol_defs["B"].modification_defs
    assert not mol_defs["B"].localization_def
    assert not mol_defs["B"].association_defs


def test_single_deg_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_deg_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    assert not mol_defs["A"].modification_defs
    assert not mol_defs["A"].localization_def
    assert not mol_defs["A"].association_defs

    assert not mol_defs["B"].modification_defs
    assert not mol_defs["B"].localization_def
    assert not mol_defs["B"].association_defs

def test_single_translation_reaction_mol_defs():
    # todo: let's implement this
    pass


def test_single_transcription_reaction_mol_defs():
    # todo: let's implement this
    pass


def test_single_localization_reaction_mol_defs():
    # todo: let's implement this
    pass


 ### RULE GENERATION ###
def test_specs_from_state_for_association_complement():
    reaction = rfs.reaction_from_string('A_ppi_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])
    molecule_def = rfr.molecule_defs_from_rxncon(rxncon_sys)

    association_specs = rfr.specs_from_state(molecule_def,rfs.state_from_string("A--B"), reaction, True)

    expected_association_spec_A = rbm.AssociationSpecification(molecule_def["A"].association_def_by_domain_name("AssocB"), rbm.OccupationStatus.not_occupied)
    expected_association_spec_B = rbm.AssociationSpecification(molecule_def["B"].association_def_by_domain_name("AssocA"), rbm.OccupationStatus.not_occupied)
    assert association_specs["A"] == expected_association_spec_A
    assert association_specs["B"] == expected_association_spec_B


def test_specs_from_state_for_association_not_complement():
    reaction = rfs.reaction_from_string('A_ppi_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])
    molecule_def = rfr.molecule_defs_from_rxncon(rxncon_sys)

    association_specs = rfr.specs_from_state(molecule_def,rfs.state_from_string("A--B"), reaction, False)

    expected_association_spec_A = rbm.AssociationSpecification(molecule_def["A"].association_def_by_domain_name("AssocB"), rbm.OccupationStatus.occupied_known_partner)
    expected_association_spec_B = rbm.AssociationSpecification(molecule_def["B"].association_def_by_domain_name("AssocA"), rbm.OccupationStatus.occupied_known_partner)
    assert association_specs["A"] == expected_association_spec_A
    assert association_specs["B"] == expected_association_spec_B


def test_specs_from_state_modification_not_complement():
    reaction = rfs.reaction_from_string('A_p+_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    molecule_def = rfr.molecule_defs_from_rxncon(rxncon_sys)

    modification_specs = rfr.specs_from_state(molecule_def,rfs.state_from_string("B-{P}"), reaction, False)

    expected_modification_spec = rbm.ModificationSpecification(molecule_def["B"].modification_def_by_domain_name("ModA"),"p")

    assert modification_specs["B"] == expected_modification_spec


def test_specs_from_state_modification_complement():
    reaction = rfs.reaction_from_string('A_p+_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    molecule_def = rfr.molecule_defs_from_rxncon(rxncon_sys)

    modification_specs = rfr.specs_from_state(molecule_def,rfs.state_from_string("B-{P}"), reaction, True)

    expected_modification_spec = rbm.ModificationSpecification(molecule_def["B"].modification_def_by_domain_name("ModA"),"u")

    assert modification_specs["B"] == expected_modification_spec


def test_specs_pair_from_flow_association():
    reaction = rfs.reaction_from_string('A_ppi_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    molecule_def = rfr.molecule_defs_from_rxncon(rxncon_sys)

    boolean_flow = flo.boolean_state_flows(reaction, rxncon_sys.strict_contingencies_for_reaction(reaction),
                                           rxncon_sys.source_contingencies_for_reaction(reaction))
    #todo: remove reaction as input
    source_specs, target_specs = rfr.specs_pair_from_flow(molecule_def, boolean_flow[0], reaction)

    expected_source_association_spec_A = {rbm.AssociationSpecification(molecule_def["A"].association_def_by_domain_name("AssocB"), rbm.OccupationStatus.not_occupied)}
    expected_source_association_spec_B = {rbm.AssociationSpecification(molecule_def["B"].association_def_by_domain_name("AssocA"), rbm.OccupationStatus.not_occupied)}

    expected_target_association_spec_A = {rbm.AssociationSpecification(molecule_def["A"].association_def_by_domain_name("AssocB"), rbm.OccupationStatus.occupied_known_partner)}
    expected_target_association_spec_B = {rbm.AssociationSpecification(molecule_def["B"].association_def_by_domain_name("AssocA"), rbm.OccupationStatus.occupied_known_partner)}

    assert source_specs == {"A": expected_source_association_spec_A, "B": expected_source_association_spec_B}
    assert target_specs == {"A": expected_target_association_spec_A, "B": expected_target_association_spec_B}

def test_reactants_from_specs_pair():
    reaction = rfs.reaction_from_string('A_p+_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    molecule_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    boolean_flow = flo.boolean_state_flows(reaction, rxncon_sys.strict_contingencies_for_reaction(reaction),
                                           rxncon_sys.source_contingencies_for_reaction(reaction))
    #todo: remove reaction as input
    source_specs, target_specs = rfr.specs_pair_from_flow(molecule_defs, boolean_flow[0], reaction)

    source_expected_reactant_A = rbm.MoleculeReactant(rbm.MoleculeSpecification(molecule_defs["A"], [], [], None))
    source_expected_reactant_B = rbm.MoleculeReactant(rbm.MoleculeSpecification(molecule_defs["B"],
                                                           [rbm.ModificationSpecification(molecule_defs["B"].modification_def_by_domain_name("ModA"),"u")],
                                                           [], None))

    target_expected_reactant_A = rbm.MoleculeReactant(rbm.MoleculeSpecification(molecule_defs["A"], [], [], None))
    target_expected_reactant_B = rbm.MoleculeReactant(rbm.MoleculeSpecification(molecule_defs["B"],
                                                           [rbm.ModificationSpecification(molecule_defs["B"].modification_def_by_domain_name("ModA"),"p")],
                                                           [], None))
    reactants = rfr.reactants_from_specs_pair(molecule_defs, (source_specs, target_specs))

    assert len(reactants[0]) == 2
    assert reactants[0][0] != reactants[0][1]
    assert reactants[0][0] in [source_expected_reactant_A, source_expected_reactant_B]
    assert reactants[0][1] in [source_expected_reactant_A, source_expected_reactant_B]

    assert len(reactants[1]) == 2
    assert reactants[1][0] in [target_expected_reactant_A, target_expected_reactant_B]
    assert reactants[1][1] in [target_expected_reactant_A, target_expected_reactant_B]
    assert reactants[1][0] != reactants[1][1]


def test_rule_from_flow_and_molecule_definitions_and_reaction():
    reaction = rfs.reaction_from_string('A_p+_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    molecule_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    boolean_flow = flo.boolean_state_flows(reaction, rxncon_sys.strict_contingencies_for_reaction(reaction),
                                           rxncon_sys.source_contingencies_for_reaction(reaction))

    left_hand_side = [rbm.MoleculeReactant(rbm.MoleculeSpecification(molecule_defs["A"], [], [], None)),
                      rbm.MoleculeReactant(rbm.MoleculeSpecification(molecule_defs["B"],
                                                                     [rbm.ModificationSpecification(molecule_defs["B"].modification_def_by_domain_name("ModA"),"u")],
                                                                     [], None))]

    right_hand_side = [rbm.MoleculeReactant(rbm.MoleculeSpecification(molecule_defs["A"], [], [], None)),
                       rbm.MoleculeReactant(rbm.MoleculeSpecification(molecule_defs["B"],
                                                                      [rbm.ModificationSpecification(molecule_defs["B"].modification_def_by_domain_name("ModA"),"p")],
                                                                      [], None))]


    rule = rfr.rule_from_flow_and_molecule_definitions_and_reaction(molecule_defs, boolean_flow[0], reaction, 1)

    assert len(rule.left_hand_side) == 2
    assert rule.left_hand_side[0] in left_hand_side
    assert rule.left_hand_side[1] in left_hand_side
    assert rule.left_hand_side[0] != rule.left_hand_side[1]

    assert len(rule.right_hand_side) == 2
    assert rule.right_hand_side[0] in right_hand_side
    assert rule.right_hand_side[1] in right_hand_side
    assert rule.right_hand_side[0] != rule.right_hand_side[1]

    assert rule.arrow_type == rbm.Arrow.irreversible
    assert rule.rates == [rbm.Parameter("k1","1")]


def test_rule_based_model_from_rxncon():
    reaction = rfs.reaction_from_string('A_ppi_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    rfr.rule_based_model_from_rxncon(rxncon_sys)



# def test_test():
#     # todo: does not work contingencies are not considered in mol_def generation
#     quick_test = qui.Quick("A_ppi_B; ! A--C")
#     molecule_def = rfr.molecule_defs_from_rxncon(quick_test.rxncon_system)
#     reaction = quick_test.rxncon_system.reactions[0]
#     boolean_flow = flo.boolean_state_flows(reaction, quick_test.rxncon_system.strict_contingencies_for_reaction(reaction),
#                                            quick_test.rxncon_system.source_contingencies_for_reaction(reaction))
#     #todo: remove reaction as input
#     source_specs, target_specs = rfr.specs_pair_from_flow(molecule_def, boolean_flow[0], quick_test.rxncon_system.reactions[0])
#
#     pass

