import rxncon.core.reaction as rxn
import rxncon.core.rxncon_system as rxs
import rxncon.core.contingency as con
import rxncon.core.effector as eff
import rxncon.syntax.rxncon_from_string as rfs
import rxncon.simulation.rule_based.rule_based_model_from_rxncon as rfr


### MOLECULE DEFINITIONS FROM SINGLE REACTION SYSTEMS ###
def test_single_ppi_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_ppi_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    mol_def_a = [x for x in mol_defs if x.name == 'A'][0]
    mol_def_b = [x for x in mol_defs if x.name == 'B'][0]

    assert not mol_def_a.modification_defs
    assert not mol_def_a.localization_def

    assert len(mol_def_a.association_defs) == 1
    assert mol_def_a.association_defs[0].domain == 'AssocB'

    assert not mol_def_b.modification_defs
    assert not mol_def_b.localization_def

    assert len(mol_def_b.association_defs) == 1
    assert mol_def_b.association_defs[0].domain == 'AssocA'


def test_single_ipi_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_[n]_ipi_A_[d]')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    mol_def_a = [x for x in mol_defs if x.name == 'A'][0]
    assert not mol_def_a.modification_defs
    assert not mol_def_a.localization_def

    assert len(mol_def_a.association_defs) == 2
    assert mol_def_a.association_defs[0] != mol_def_a.association_defs[1]
    assert mol_def_a.association_defs[0].domain in ["n", "d"]
    assert mol_def_a.association_defs[1].domain in ["n", "d"]


def test_single_modification_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_p+_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    mol_def_a = [x for x in mol_defs if x.name == 'A'][0]
    mol_def_b = [x for x in mol_defs if x.name == 'B'][0]

    assert not mol_def_a.association_defs
    assert not mol_def_a.modification_defs
    assert not mol_def_a.localization_def

    assert not mol_def_b.association_defs
    assert not mol_def_b.localization_def

    assert len(mol_def_b.modification_defs) == 1
    assert mol_def_b.modification_defs[0].domain == "ModA"
    assert mol_def_b.modification_defs[0].valid_modifiers == ["u","p"]


def test_single_synth_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_syn_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    mol_def_a = [x for x in mol_defs if x.name == 'A'][0]
    mol_def_b = [x for x in mol_defs if x.name == 'B'][0]

    assert not mol_def_a.modification_defs
    assert not mol_def_a.localization_def
    assert not mol_def_a.association_defs

    assert not mol_def_b.modification_defs
    assert not mol_def_b.localization_def
    assert not mol_def_b.association_defs


def test_single_deg_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_deg_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    mol_def_a = [x for x in mol_defs if x.name == 'A'][0]
    mol_def_b = [x for x in mol_defs if x.name == 'B'][0]

    assert not mol_def_a.modification_defs
    assert not mol_def_a.localization_def
    assert not mol_def_a.association_defs

    assert not mol_def_b.modification_defs
    assert not mol_def_b.localization_def
    assert not mol_def_b.association_defs

def test_single_translation_reaction_mol_defs():
    # todo: let's implement this
    pass


def test_single_transcription_reaction_mol_defs():
    # todo: let's implement this
    pass


def test_single_localization_reaction_mol_defs():
    # todo: let's implement this
    pass



### SET STUFF ###
def test_set_from_contingencies():
    reaction = rfs.reaction_from_string('A_ppi_B')
    state1 = rfs.state_from_string('A_[x]-{p}')
    state2 = rfs.state_from_string('B_[y]-{p}')
    state3 = rfs.state_from_string('A_[z]-{p}')
    state4 = rfs.state_from_string('B_[w]-{p}')

    contingency = con.Contingency(reaction,
                                  con.ContingencyType.requirement,
                                  eff.AndEffector(eff.OrEffector(eff.StateEffector(state1),
                                                                 eff.StateEffector(state2)),
                                                  eff.OrEffector(eff.StateEffector(state3),
                                                                 eff.StateEffector(state4))))

    source_contingency = con.Contingency(reaction,
                                         con.ContingencyType.inhibition,
                                         eff.StateEffector(rfs.state_from_string('A--B')))

    the_set = rfr.set_from_contingencies([contingency])

    rule_conditions = rfr.base_rule_conditions_from_strict_and_source_contingencies([contingency], [source_contingency])

    for rule in rule_conditions:
        print(rule)










