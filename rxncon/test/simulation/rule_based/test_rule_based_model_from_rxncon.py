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

def test_double_ppi_reaction_mol_defs():
    reaction1 = rfs.reaction_from_string('A_ppi_B')
    reaction2 = rfs.reaction_from_string('A_ppi_C')
    rxncon_sys = rxs.RxnConSystem([reaction1, reaction2], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    mol_def_a = [x for x in mol_defs if x.name == 'A'][0]
    mol_def_b = [x for x in mol_defs if x.name == 'B'][0]
    mol_def_c = [x for x in mol_defs if x.name == 'C'][0]

    assert not mol_def_a.modification_defs
    assert not mol_def_a.localization_def

    assert len(mol_def_a.association_defs) == 2
    assert mol_def_a.association_defs[0].domain != mol_def_a.association_defs[1].domain
    assert mol_def_a.association_defs[0].domain in ['AssocB', 'AssocC']
    assert mol_def_a.association_defs[1].domain in ['AssocB', 'AssocC']

    assert not mol_def_b.modification_defs
    assert not mol_def_b.localization_def

    assert len(mol_def_b.association_defs) == 1
    assert mol_def_b.association_defs[0].domain == 'AssocA'

    assert not mol_def_c.modification_defs
    assert not mol_def_c.localization_def

    assert len(mol_def_c.association_defs) == 1
    assert mol_def_c.association_defs[0].domain == 'AssocA'


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

def test_double_modification_reaction_mol_defs():
    reaction = rfs.reaction_from_string('A_p+_B')
    reaction1 = rfs.reaction_from_string('C_p+_B_[c]')
    rxncon_sys = rxs.RxnConSystem([reaction, reaction1], [])

    mol_defs = rfr.molecule_defs_from_rxncon(rxncon_sys)

    mol_def_a = [x for x in mol_defs if x.name == 'A'][0]
    mol_def_b = [x for x in mol_defs if x.name == 'B'][0]
    mol_def_c = [x for x in mol_defs if x.name == 'C'][0]

    assert not mol_def_a.association_defs
    assert not mol_def_a.modification_defs
    assert not mol_def_a.localization_def

    assert not mol_def_c.association_defs
    assert not mol_def_c.modification_defs
    assert not mol_def_c.localization_def

    assert not mol_def_b.association_defs
    assert not mol_def_b.localization_def

    assert len(mol_def_b.modification_defs) == 2
    assert mol_def_b.modification_defs[0].domain != mol_def_b.modification_defs[1].domain
    assert mol_def_b.modification_defs[0].domain in ["ModA","c"]
    assert mol_def_b.modification_defs[1].domain in ["ModA","c"]

    assert mol_def_b.modification_defs[0].valid_modifiers == ["u","p"]
    assert mol_def_b.modification_defs[1].valid_modifiers == ["u","p"]



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


 ### RULE GENERATION ###

def test_rule_based_model_from_rxncon():
    reaction = rfs.reaction_from_string('A_ppi_B')
    rxncon_sys = rxs.RxnConSystem([reaction], [])

    #rfr.rule_based_model_from_rxncon(rxncon_sys)




