import os

import rxncon.input.sbtab.sbtab as sbt


TIGER_FILENAME = 'Tiger_et_al_TableS1_SBtab_ContingencyID.tsv'
TIGER_PATH = os.path.join(os.path.dirname(__file__), TIGER_FILENAME)


def test_sbtab_from_list_of_lists():
    sbtab_input = [
        ['!!SBtab SBtabVersion "0.8" Document="Hynne_2001" TableType="Quantity" TableName="Quantity-parameters"'],
        ['!Quantity', '!Name', '!Compound', '!Reaction', '!Location', '!Value', '!Unit', '!SBOTerm', '!QuantityType'],
        ['Par1', 'k1', '', 'vGlcTrans', 'plasmamembrane', '0.14', '1/s', 'SBO:0000022', 'forward unimolecular rate constant'],
        ['Par2', 'Km', 'Glc', 'vHK', 'cytosol', '0.31', 'mM', 'SBO:0000322', 'Michaelis constant for substrate'],
        ['Par3', 'Vm', '', 'vHK', 'cytosol', '2.14', 'mM/s', 'SBO:0000153', 'forward rate constant']
    ]

    sbtab = sbt.SBtabData(sbtab_input)

    assert len(sbtab.entries) == 3

    entry = sbtab.entries[0]
    assert entry.Quantity == 'Par1'
    assert entry.Name == 'k1'
    assert entry.Compound == ''
    assert entry.Reaction == 'vGlcTrans'
    assert entry.Location == 'plasmamembrane'
    assert entry.Value == '0.14'
    assert entry.Unit == '1/s'
    assert entry.SBOTerm == 'SBO:0000022'
    assert entry.QuantityType == 'forward unimolecular rate constant'


def test_sbtab_tiger_contingencies():
    sbtab = sbt.sbtab_data_from_file(TIGER_PATH)

    assert len(sbtab.entries) == 313
    print(dir(sbtab.entries[0]))
