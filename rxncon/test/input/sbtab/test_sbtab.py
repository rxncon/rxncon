
import rxncon.input.sbtab.sbtab as sbt


def test_sbtab_from_list_of_lists():
    sbtab_input = [
        ['!!SBtab SBtabVersion "0.8" Document="Hynne_2001" TableType="Quantity" TableName="Quantity-parameters"'],
        ['!Quantity', '!Name', '!Compound', '!Reaction', '!Location', '!Value', '!Unit', '!SBOTerm', '!QuantityType'],
        ['Par1', 'k1', '', 'vGlcTrans', 'plasmamembrane', '0.14', '1/s', 'SBO:0000022', 'forward unimolecular rate constant'],
        ['Par2', 'Km', 'Glc', 'vHK', 'cytosol', '0.31', 'mM', 'SBO:0000322', 'Michaelis constant for substrate'],
        ['Par3', 'Vm', '', 'vHK', 'cytosol', '2.14', 'mM/s', 'SBO:0000153', 'forward rate constant']
    ]

    sbtan = sbt.SBtabData(sbtab_input)