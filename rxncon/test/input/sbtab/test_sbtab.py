import os
import math
import pytest

from rxncon.input.sbtab.sbtab import sbtab_data_from_file, ValidatedSBtabData, SBtabData


TIGER_CONTINGENCY_FILENAME = 'Tiger_et_al_TableS1_SBtab_ContingencyID.tsv'
TIGER_CONTINGENCY_PATH     = os.path.join(os.path.dirname(__file__), TIGER_CONTINGENCY_FILENAME)

DEFINITIONS_FILENAME = 'definitions.tsv'
DEFINITIONS_PATH     = os.path.join(os.path.dirname(__file__), DEFINITIONS_FILENAME)

RXNCON_DEFINITIONS_FILENAME = 'rxncon_Definition.tsv'
RXNCON_DEFINITIONS_PATH     = os.path.join(os.path.dirname(__file__), RXNCON_DEFINITIONS_FILENAME)


def test_sbtab_definition_from_file() -> None:
    sbtab = sbtab_data_from_file(DEFINITIONS_PATH)
    assert len(sbtab.entries) == 249


def test_sbtab_header_version_ambiguity_parses() -> None:
    input_version_separator_space = [
        ['!!SBtab SBtabVersion "0.8" Document="Hynne_2001" TableType="Quantity" TableName="Quantity-parameters"'],
        ['!Quantity', '!Name', '!Compound', '!Reaction', '!Location', '!Value', '!Unit', '!SBOTerm', '!QuantityType'],
        ['Par1', 'k1', '', 'vGlcTrans', 'plasmamembrane', '0.14', '1/s', 'SBO:0000022', 'forward unimolecular rate constant']
    ]

    input_version_separator_equals = [
        ['!!SBtab SBtabVersion="0.8" Document="Hynne_2001" TableType="Quantity" TableName="Quantity-parameters"'],
        ['!Quantity', '!Name', '!Compound', '!Reaction', '!Location', '!Value', '!Unit', '!SBOTerm', '!QuantityType'],
        ['Par1', 'k1', '', 'vGlcTrans', 'plasmamembrane', '0.14', '1/s', 'SBO:0000022', 'forward unimolecular rate constant']
    ]

    sbtab_space  = SBtabData(input_version_separator_space)
    sbtab_equals = SBtabData(input_version_separator_equals)

    assert sbtab_space.version == sbtab_equals.version == '0.8'


def test_validated_sbtab_from_list_of_lists() -> None:
    sbtab_input = [
        ['!!SBtab SBtabVersion "0.8" Document="Hynne_2001" TableType="Quantity" TableName="Quantity-parameters"'],
        ['!Quantity', '!Name', '!Compound', '!Reaction', '!Location', '!Value', '!Unit', '!SBOTerm', '!QuantityType'],
        ['Par1', 'k1', '', 'vGlcTrans', 'plasmamembrane', '0.14', '1/s', 'SBO:0000022', 'forward unimolecular rate constant'],
        ['Par2', 'Km', 'Glc', 'vHK', 'cytosol', '0.31', 'mM', 'SBO:0000322', 'Michaelis constant for substrate'],
        ['Par3', 'Vm', '', 'vHK', 'cytosol', '2.14', 'mM/s', 'SBO:0000153', 'forward rate constant']
    ]

    definitions = sbtab_data_from_file(DEFINITIONS_PATH)
    sbtab = ValidatedSBtabData(sbtab_input, definitions)

    assert sbtab.version       == '0.8'
    assert sbtab.document_name == 'Hynne_2001'
    assert sbtab.table_type    == 'Quantity'
    assert sbtab.table_name    == 'Quantity-parameters'
    assert len(sbtab.entries)  == 3

    # Assert both "entry[property]" and the "entry.property" work.
    entry = sbtab.entries[0]
    assert entry.field_names == ['Quantity', 'Name', 'Compound', 'Reaction', 'Location', 'Value', 'Unit', 'SBOTerm', 'QuantityType']  # type: ignore
    assert [entry[field] for field in entry.field_names if field != 'Value'] == ['Par1', 'k1', '', 'vGlcTrans', 'plasmamembrane', '1/s', 'SBO:0000022', 'forward unimolecular rate constant']  # type: ignore

    assert entry.Quantity     == 'Par1'                                 # type: ignore
    assert entry.Name         == 'k1'                                   # type: ignore
    assert entry.Compound     == ''                                     # type: ignore
    assert entry.Reaction     == 'vGlcTrans'                            # type: ignore
    assert entry.Location     == 'plasmamembrane'                       # type: ignore
    assert math.isclose(entry.Value, 0.14)                              # type: ignore
    assert entry.Unit         == '1/s'                                  # type: ignore
    assert entry.SBOTerm      == 'SBO:0000022'                          # type: ignore
    assert entry.QuantityType == 'forward unimolecular rate constant'   # type: ignore


def test_validated_sbtab_raises_value_error_wrong_type() -> None:
    sbtab_input = [
        ['!!SBtab SBtabVersion "0.8" Document="Hynne_2001" TableType="Quantity" TableName="Quantity-parameters"'],
        ['!Quantity', '!Name', '!Compound', '!Reaction', '!Location', '!Value', '!Unit', '!SBOTerm', '!QuantityType'],
        ['Par1', 'k1', '', 'vGlcTrans', 'plasmamembrane', 'this_should_be_a_float', '1/s', 'SBO:0000022', 'forward unimolecular rate constant']
    ]

    definitions = sbtab_data_from_file(DEFINITIONS_PATH)

    with pytest.raises(ValueError):
        ValidatedSBtabData(sbtab_input, definitions)


def test_input_type_validated_sbtab_() -> None:

    sbtab_input= [
        ['!!SBtab SBtabVersion "0.8" Document="Hynne_2001" TableType="Compound" TableName="Test_Compound"'],
        ['!Charge', '!Compound', '!InitialValue', '!IsConstant', '!Comment'],
        ['1',  'compound', '1.5', 'True', 'some mode comments']
    ]

    definitions = sbtab_data_from_file(DEFINITIONS_PATH)

    sbtab = ValidatedSBtabData(sbtab_input, definitions)

    assert sbtab.version       == '0.8'
    assert sbtab.document_name == 'Hynne_2001'
    assert sbtab.table_type    == 'Compound'
    assert sbtab.table_name    == 'Test_Compound'
    assert len(sbtab.entries)  == 1

    # Assert both "entry[property]" and the "entry.property" work.
    entry = sbtab.entries[0]
    assert entry.field_names   == ['Charge', 'Compound', 'InitialValue', 'IsConstant', 'Comment']  # type: ignore

    assert entry.Charge        == 1                                   # type: ignore
    assert entry.Compound      == 'compound'                               # type: ignore
    assert math.isclose(entry.InitialValue, 1.5)                        # type: ignore
    assert entry.IsConstant    == True                          # type: ignore
    assert entry.Comment       == 'some mode comments'                        # type: ignore


def test_validated_sbtab_tiger_contingencies() -> None:
    definitions = sbtab_data_from_file(RXNCON_DEFINITIONS_PATH)

    assert len(definitions.entries) > 0
    assert definitions.entries[0].__class__.__name__ == 'RxnconDefinitionsTable'

    sbtab = sbtab_data_from_file(TIGER_CONTINGENCY_PATH, definitions=definitions)
    assert len(sbtab.entries) == 313

    entry = sbtab.entries[0]
    assert entry.ContingencyID     == 'ID1'                                                     # type: ignore
    assert entry.Target            == 'Hog1_[KD]_P+_Rck2_[c(S519)]'                             # type: ignore
    assert entry.Contingency       == '!'                                                       # type: ignore
    assert entry.Modifier          == 'Hog1_[(T174)]-{P}'                                       # type: ignore
    assert entry.PubMedIdentifiers == '10805732'                                                # type: ignore
    assert entry.Quality           == 'in vitro kinase assay'                                   # type: ignore
    assert entry.Comment           == 'Hog1 activation by constitutive active Pbs2 needed.'     # type: ignore
