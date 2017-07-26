from rxncon.input.quick.quick import Quick
from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.simulation.boolean.boolean_model import boolean_model_from_rxncon, ReactionTarget, \
    StateTarget, SmoothingStrategy
from rxncon.simulation.boolean.boolnet_from_boolean_model import boolnet_from_boolean_model
from rxncon.test.simulation.boolean.utils import target_from_str
from rxncon.venntastic.sets import venn_from_str


def test_simple_system() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                                    A_[b]_ppi-_B_[a]
                                                    C_p+_A_[(r)]
                                                    D_p-_A_[(r)]""").rxncon_system)

    # Component expressions.
    A = '(( A_[b]--0 | A_[b]--B_[a] ) & ( A_[(r)]-{0} | A_[(r)]-{p} ))'
    B = '( B_[a]--0 | A_[b]--B_[a] )'
    C = 'C'
    D = 'D'

    expected_rules = {
        'A_[b]_ppi+_B_[a]': '{0} & {1} & A_[(r)]-{{p}}'.format(A, B),
        'A_[b]_ppi-_B_[a]': '{0} & {1}'.format(A, B),
        'C_p+_A_[(r)]': '{0} & {1}'.format(A, C),
        'D_p-_A_[(r)]': '{0} & {1}'.format(A, D),
        'A_[(r)]-{p}': '{0} & (( C_p+_A_[(r)] & A_[(r)]-{{0}} ) | '
                       '( A_[(r)]-{{p}} & ~( D_p-_A_[(r)] & A_[(r)]-{{p}} )))'.format(A),
        'A_[(r)]-{0}': '{0} & (( D_p-_A_[(r)] & A_[(r)]-{{p}} ) | '
                       '( A_[(r)]-{{0}} & ~( C_p+_A_[(r)] & A_[(r)]-{{0}} )))'.format(A),
        'A_[b]--0': '{0} & (( A_[b]_ppi-_B_[a] & A_[b]--B_[a] ) | '
                    '( A_[b]--0 & ~( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 )))'.format(A),
        'A_[b]--B_[a]': '{0} & {1} & (( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 ) | '
                        '( A_[b]--B_[a] & ~( A_[b]_ppi-_B_[a] & A_[b]--B_[a] )))'.format(A, B),
        'B_[a]--0': '{0} & (( A_[b]_ppi-_B_[a] & A_[b]--B_[a] ) | '
                    '( B_[a]--0 & ~( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 )))'.format(B),
        'C': '{0}'.format(C),
        'D': '{0}'.format(D),
    }

    assert len(boolean_model.update_rules) == len(expected_rules)

    for update_rule in boolean_model.update_rules:
        assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                 target_from_str))


def test_simple_system_with_input_state() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p} ; ! [BLAAT]
                                                        A_[b]_ppi-_B_[a]
                                                        C_p+_A_[(r)]
                                                        D_p-_A_[(r)]""").rxncon_system)

    # Component expressions.
    A = '(( A_[b]--0 | A_[b]--B_[a] ) & ( A_[(r)]-{0} | A_[(r)]-{p} ))'
    B = '( B_[a]--0 | A_[b]--B_[a] )'
    C = 'C'
    D = 'D'

    expected_rules = {
        'A_[b]_ppi+_B_[a]': '{0} & {1} & A_[(r)]-{{p}} & [BLAAT]'.format(A, B),
        'A_[b]_ppi-_B_[a]': '{0} & {1}'.format(A, B),
        'C_p+_A_[(r)]': '{0} & {1}'.format(A, C),
        'D_p-_A_[(r)]': '{0} & {1}'.format(A, D),
        'A_[(r)]-{p}': '{0} & (( C_p+_A_[(r)] & A_[(r)]-{{0}} ) | '
                       '( A_[(r)]-{{p}} & ~( D_p-_A_[(r)] & A_[(r)]-{{p}} )))'.format(A),
        'A_[(r)]-{0}': '{0} & (( D_p-_A_[(r)] & A_[(r)]-{{p}} ) | '
                       '( A_[(r)]-{{0}} & ~( C_p+_A_[(r)] & A_[(r)]-{{0}} )))'.format(A),
        'A_[b]--0': '{0} & (( A_[b]_ppi-_B_[a] & A_[b]--B_[a] ) | '
                    '( A_[b]--0 & ~( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 )))'.format(A),
        'A_[b]--B_[a]': '{0} & {1} & (( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 ) | '
                        '( A_[b]--B_[a] & ~( A_[b]_ppi-_B_[a] & A_[b]--B_[a] )))'.format(A, B),
        'B_[a]--0': '{0} & (( A_[b]_ppi-_B_[a] & A_[b]--B_[a] ) | '
                    '( B_[a]--0 & ~( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 )))'.format(B),
        'C': '{0}'.format(C),
        'D': '{0}'.format(D),
        '[BLAAT]': '[BLAAT]'
    }

    assert len(boolean_model.update_rules) == len(expected_rules)

    for update_rule in boolean_model.update_rules:
        assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                 target_from_str))


def test_trsl_trsc_deg() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""PolII_trsc_TargetGene
                                                       Ribo_trsl_TargetmRNA
                                                       Nuclease_deg_TargetmRNA
                                                       Protease_deg_Target""").rxncon_system)

    expected_rules = {
        'PolII_trsc_TargetGene': '( PolII & TargetGene )',
        'Ribo_trsl_TargetmRNA': '( Ribo & TargetmRNA )',
        'Nuclease_deg_TargetmRNA': '( Nuclease & TargetmRNA )',
        'Protease_deg_Target': '( Protease & Target )',
        'PolII': '( PolII )',
        'TargetGene': '( TargetGene )',
        'Ribo': '( Ribo )',
        'TargetmRNA': '( PolII_trsc_TargetGene | ( TargetmRNA & ~( Nuclease_deg_TargetmRNA )))',
        'Nuclease': '( Nuclease )',
        'Protease': '( Protease )',
        'Target': '( Ribo_trsl_TargetmRNA | ( Target & ~( Protease_deg_Target )))'
    }

    assert len(boolean_model.update_rules) == len(expected_rules)

    for update_rule in boolean_model.update_rules:
        assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                 target_from_str))


def test_deg_without_contingency() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""B_p+_A_[(res)]
                                                       D_ub+_A_[(res)]
                                                       D_p+_A_[(other)]
                                                       C_deg_A""").rxncon_system)

    # Component expression.
    A = '(( A_[(res)]-{0} | A_[(res)]-{p} | A_[(res)]-{ub} ) & ( A_[(other)]-{0} | A_[(other)]-{p} ))'

    expected_rules = {
        'B_p+_A_[(res)]': '( B & {} )'.format(A),
        'D_ub+_A_[(res)]': '( D & {} )'.format(A),
        'D_p+_A_[(other)]': '( D & {} )'.format(A),
        'C_deg_A': '( C & {} )'.format(A),
        'A_[(res)]-{0}': '( {} & ~( C_deg_A ) & A_[(res)]-{{0}} & ~( B_p+_A_[(res)] & A_[(res)]-{{0}} ) & '
                         '~( D_ub+_A_[(res)] & A_[(res)]-{{0}} ))'.format(A),
        'A_[(res)]-{p}': '( {} & ~( C_deg_A ) & (( B_p+_A_[(res)] & A_[(res)]-{{0}} ) | A_[(res)]-{{p}} ))'.format(A),
        'A_[(res)]-{ub}': '( {} & ~( C_deg_A ) & (( D_ub+_A_[(res)] & A_[(res)]-{{0}} ) | A_[(res)]-{{ub}} ))'.format(
            A),
        'A_[(other)]-{0}': '( {} & ~( C_deg_A ) & A_[(other)]-{{0}} & ~( D_p+_A_[(other)] & '
                           'A_[(other)]-{{0}} ))'.format(A),
        'A_[(other)]-{p}': '( {} & ~( C_deg_A ) & (( D_p+_A_[(other)] & A_[(other)]-{{0}} ) | '
                           'A_[(other)]-{{p}} ))'.format(A),
        'B': 'B',
        'C': 'C',
        'D': 'D'
    }

    assert len(boolean_model.update_rules) == len(expected_rules)

    for update_rule in boolean_model.update_rules:
        assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                 target_from_str))


def test_deg_with_requirement() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""B_p+_A_[(res)]
                                                       D_ub+_A_[(res)]
                                                       D_p+_A_[(other)]
                                                       C_deg_A; ! A_[(res)]-{p}""").rxncon_system)

    # Component expression.
    A = '(( A_[(res)]-{0} | A_[(res)]-{p} | A_[(res)]-{ub} ) & ( A_[(other)]-{0} | A_[(other)]-{p} ))'

    expected_rules = {
        'B_p+_A_[(res)]': '( B & {} )'.format(A),
        'D_ub+_A_[(res)]': '( D & {} )'.format(A),
        'D_p+_A_[(other)]': '( D & {} )'.format(A),
        'C_deg_A': '( C & {} & A_[(res)]-{{p}} )'.format(A),
        'A_[(res)]-{0}': '( {} & A_[(res)]-{{0}} & ~( B_p+_A_[(res)] & A_[(res)]-{{0}} ) & '
                         '~( D_ub+_A_[(res)] & A_[(res)]-{{0}} ))'.format(A),
        'A_[(res)]-{p}': '( {} & ~( C_deg_A ) & (( B_p+_A_[(res)] & A_[(res)]-{{0}} ) | A_[(res)]-{{p}} ))'.format(A),
        'A_[(res)]-{ub}': '( {} & (( D_ub+_A_[(res)] & A_[(res)]-{{0}} ) | A_[(res)]-{{ub}} ))'.format(A),
        'A_[(other)]-{0}': '( {} & A_[(other)]-{{0}} & ~( D_p+_A_[(other)] & A_[(other)]-{{0}} ))'.format(A),
        'A_[(other)]-{p}': '( {} & (( D_p+_A_[(other)] & A_[(other)]-{{0}} ) | A_[(other)]-{{p}} ))'.format(A),
        'B': 'B',
        'C': 'C',
        'D': 'D'
    }

    assert len(boolean_model.update_rules) == len(expected_rules)

    for update_rule in boolean_model.update_rules:
        assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                 target_from_str))


def test_deg_with_contingency_on_subject() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_deg_B ; ! A-{p}
                                                    C_p+_A""").rxncon_system)

    assert len(boolean_model.update_rules) == 6

    found_B, found_C = False, False

    for update_rule in boolean_model.update_rules:
        if update_rule.target == target_from_str('B'):
            found_B = True
            assert update_rule.factor.is_equivalent_to(venn_from_str('B & ~( A_deg_B )', target_from_str))
        elif update_rule.target == target_from_str('C'):
            found_C = True
            assert update_rule.factor.is_equivalent_to(venn_from_str('C', target_from_str))

    assert found_B and found_C


def test_deg_with_inhibition() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""B_p+_A_[(res)]
                                                       D_ub+_A_[(res)]
                                                       D_p+_A_[(other)]
                                                       C_deg_A; x A_[(res)]-{p}""").rxncon_system)

    assert len(boolean_model.reaction_target_by_name('C_deg_A').degraded_targets) == 2
    assert target_from_str('A_[(res)]-{0}') in boolean_model.reaction_target_by_name('C_deg_A').degraded_targets
    assert target_from_str('A_[(res)]-{ub}') in boolean_model.reaction_target_by_name('C_deg_A').degraded_targets


def test_deg_with_boolean_OR() -> None:
    """The C_deg_A reaction will be split into two ReactionTargets, one degrading A_[(r1)]-{p}, the other
    degrading A_[(r2)]-{p}. This choice which variant degrades which A is not deterministic, so the test splits
    into two branches here."""
    boolean_model = boolean_model_from_rxncon(Quick("""B_p+_A_[(r1)]
                                                       D_p+_A_[(r2)]
                                                       C_deg_A; ! <X>
                                                       <X>; OR A_[(r1)]-{p}
                                                       <X>; OR A_[(r2)]-{p}""").rxncon_system)

    A = '(( A_[(r1)]-{0} | A_[(r1)]-{p} ) & ( A_[(r2)]-{0} | A_[(r2)]-{p} ))'

    target_to_factor = {rule.target: rule.factor for rule in boolean_model.update_rules}

    # C_deg_A#c0 degrades A_[(r1)]-{p}, C_deg_A#c1 degrades A_[(r2)]-{p}
    if target_to_factor[target_from_str('C_deg_A#c0')].is_equivalent_to(
            venn_from_str('( {} & C & A_[(r1)]-{{p}} )'.format(A),
                          target_from_str)):
        expected_rules = {
            'C_deg_A#c0': '( {} & C & A_[(r1)]-{{p}} )'.format(A),
            'C_deg_A#c1': '( {} & C & A_[(r2)]-{{p}} )'.format(A),
            'A_[(r1)]-{p}': '( {} & ~( C_deg_A#c0 ) & ( A_[(r1)]-{{p}} | ( B_p+_A_[(r1)] & A_[(r1)]-{{0}} )))'.format(
                A),
            'A_[(r2)]-{p}': '( {} & ~( C_deg_A#c1 ) & ( A_[(r2)]-{{p}} | ( D_p+_A_[(r2)] & A_[(r2)]-{{0}} )))'.format(
                A),
            'B': 'B',
            'C': 'C',
            'D': 'D'
        }

        for target_str, factor_str in expected_rules.items():
            assert target_to_factor[target_from_str(target_str)].is_equivalent_to(venn_from_str(factor_str,
                                                                                                target_from_str))
    # C_deg_A#c0 degrades A_[(r2)]-{p}, C_deg_A#c1 degrades A_[(r1)]-{p}
    elif target_to_factor[target_from_str('C_deg_A#c0')].is_equivalent_to(
            venn_from_str('( {} & C & A_[(r2)]-{{p}} )'.format(A),
                          target_from_str)):
        expected_rules = {
            'C_deg_A#c0': '( {} & C & A_[(r2)]-{{p}} )'.format(A),
            'C_deg_A#c1': '( {} & C & A_[(r1)]-{{p}} )'.format(A),
            'A_[(r1)]-{p}': '( {} & ~( C_deg_A#c1 ) & ( A_[(r1)]-{{p}} | ( B_p+_A_[(r1)] & A_[(r1)]-{{0}} )))'.format(
                A),
            'A_[(r2)]-{p}': '( {} & ~( C_deg_A#c0 ) & ( A_[(r2)]-{{p}} | ( D_p+_A_[(r2)] & A_[(r2)]-{{0}} )))'.format(
                A),
            'B': 'B',
            'C': 'C',
            'D': 'D'
        }

        for target_str, factor_str in expected_rules.items():
            assert target_to_factor[target_from_str(target_str)].is_equivalent_to(venn_from_str(factor_str,
                                                                                                target_from_str))
    else:
        raise AssertionError


def test_deg_with_boolean_AND() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_[x]_ppi+_B_[y]
                                                       D_p+_A_[(r)]
                                                       C_deg_A; ! <x>
                                                       <x>; AND A_[x]--B_[y]; AND A_[(r)]-{p}""").rxncon_system)
    target_to_factor = {rule.target: rule.factor for rule in boolean_model.update_rules}

    expected_rules = {
        'C_deg_A': '(( C & A_[x]--B_[y] & A_[(r)]-{p} & ( A_[x]--0 | A_[x]--B_[y] ) & ( A_[(r)]-{0} | A_[(r)]-{p} )))',
        'D_p+_A_[(r)]': '(( D & ( A_[x]--0 | A_[x]--B_[y] ) & ( A_[(r)]-{0} | A_[(r)]-{p} )))',
        'A_[x]_ppi+_B_[y]': '((( A_[x]--0 | A_[x]--B_[y] ) & ( A_[(r)]-{0} | A_[(r)]-{p} ) & '
                            '( A_[x]--B_[y] | B_[y]--0 )))',
        'B_[y]--0': '((( A_[x]--B_[y] & C_deg_A & ( A_[x]--B_[y] | B_[y]--0 )) | '
                    '( B_[y]--0 & ( A_[x]--B_[y] | B_[y]--0 ) & ~(( B_[y]--0 & A_[x]_ppi+_B_[y] & A_[x]--0 )))))',
        'A_[(r)]-{p}': '((( A_[(r)]-{0} & D_p+_A_[(r)] & ( A_[x]--0 | A_[x]--B_[y] ) & '
                       '( A_[(r)]-{0} | A_[(r)]-{p} ) & ~( C_deg_A  )) | ( A_[(r)]-{p} & '
                       '( A_[x]--0 | A_[x]--B_[y] ) & ( A_[(r)]-{0} | A_[(r)]-{p} ) & ~( C_deg_A ))))',
        'A_[x]--B_[y]': '(( A_[x]--0 & B_[y]--0 & A_[x]_ppi+_B_[y] & ( A_[(r)]-{0} | '
                        'A_[(r)]-{p} ) & ( A_[x]--B_[y] | A_[x]--0 ) & ( A_[x]--B_[y] | B_[y]--0 ) & ~( C_deg_A )) | '
                        '( A_[x]--B_[y] & ( A_[(r)]-{0} | A_[(r)]-{p} ) & ( A_[x]--B_[y] | A_[x]--0 ) & ( A_[x]--B_[y] '
                        '| B_[y]--0 ) & ~(( A_[x]--B_[y] & C_deg_A )) & ~( C_deg_A )))',
        'A_[(r)]-{0}': '( A_[(r)]-{0} & ( A_[x]--0 | A_[x]--B_[y] ) & '
                       '( A_[(r)]-{0} | A_[(r)]-{p} ) & ~(( A_[(r)]-{0} & D_p+_A_[(r)] )))',
        'A_[x]--0': '(( A_[x]--0 & ( A_[x]--0 | A_[x]--B_[y] ) & '
                    '( A_[(r)]-{0} | A_[(r)]-{p} ) & ~(( A_[x]--0 & A_[x]_ppi+_B_[y] & B_[y]--0 ))))',
        'D': 'D',
        'C': 'C'
    }

    for target_str, factor_str in expected_rules.items():
        assert target_to_factor[target_from_str(target_str)].is_equivalent_to(venn_from_str(factor_str,
                                                                                            target_from_str))


def test_deg_with_boolean_NOT() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_[x]_ppi+_B_[y]
                                                       D_p+_A_[(r)]
                                                       C_deg_A; ! <x>
                                                       <x>; AND A_[x]--B_[y]; AND <NOT>
                                                       <NOT>; NOT A_[(r)]-{p}""").rxncon_system)
    target_to_factor = {rule.target: rule.factor for rule in boolean_model.update_rules}

    expected_rules = {
        'C_deg_A': '( C & A_[x]--B_[y] & ~( A_[(r)]-{p} ) & '
                   '( A_[x]--0 | A_[x]--B_[y] ) & ( A_[(r)]-{0} | A_[(r)]-{p} ))',
        'B_[y]--0': '(( A_[x]--B_[y] & C_deg_A & ( B_[y]--0 | A_[x]--B_[y] )) | '
                    '( B_[y]--0 & ( B_[y]--0 | A_[x]--B_[y] ) & ~(( B_[y]--0 & A_[x]_ppi+_B_[y] & A_[x]--0 ))))',
        'A_[(r)]-{p}': '(( A_[(r)]-{0} & ~( C_deg_A ) & D_p+_A_[(r)] & ( A_[x]--0 | A_[x]--B_[y] ) & ( A_[(r)]-{0} | '
                       'A_[(r)]-{p} )) | ( A_[(r)]-{p} & ( A_[x]--0 | A_[x]--B_[y] ) & ( A_[(r)]-{0} | A_[(r)]-{p} )))',
    }

    for target_str, factor_str in expected_rules.items():
        assert target_to_factor[target_from_str(target_str)].is_equivalent_to(
            venn_from_str(factor_str, target_from_str))


def test_deg_with_interaction() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_[x]_ppi+_B_[y]
                                                       D_p+_A_[(r)]
                                                       C_deg_A""").rxncon_system)
    target_to_factor = {rule.target: rule.factor for rule in boolean_model.update_rules}

    A = '( ( A_[x]--B_[y] | A_[x]--0 ) & ( A_[(r)]-{p} | A_[(r)]-{0} ) )'

    expected_rules = {
        'A_[x]--0': '( A_[x]--0 & ( A_[x]--B_[y] | A_[x]--0 ) & ( A_[(r)]-{p} | A_[(r)]-{0} ) & '
                    '~(( A_[x]_ppi+_B_[y] & A_[x]--0 & B_[y]--0 )) & ~( C_deg_A ))',
        'B_[y]--0': '(( C_deg_A & A_[x]--B_[y] & ( B_[y]--0 | A_[x]--B_[y] )) | ( B_[y]--0 & '
                    '( B_[y]--0 | A_[x]--B_[y] ) & ~(( A_[x]--0 & B_[y]--0 & A_[x]_ppi+_B_[y] & ~( C_deg_A )))))',
        'A_[x]--B_[y]': '{} & ~( C_deg_A ) & (( A_[x]_ppi+_B_[y] & A_[x]--0 & B_[y]--0 ) | '
                        '( A_[x]--B_[y] & ~( C_deg_A & A_[x]--B_[y] ) ) )'.format(A)
    }

    for target_str, factor_str in expected_rules.items():
        assert target_to_factor[target_from_str(target_str)].is_equivalent_to(venn_from_str(factor_str,
                                                                                            target_from_str))


def test_deg_with_interaction_as_requirement() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_[x]_ppi+_B_[y]
                                                       D_p+_A_[(r)]
                                                       C_deg_A; ! A_[x]--B_[y]""").rxncon_system)

    assert boolean_model.reaction_target_by_name('C_deg_A').degraded_targets == [target_from_str('A_[x]--B_[y]')]
    assert boolean_model.reaction_target_by_name('C_deg_A').consumed_targets == [target_from_str('A_[x]--B_[y]')]
    assert boolean_model.reaction_target_by_name('C_deg_A').produced_targets == [target_from_str('B_[y]--0')]
    assert boolean_model.reaction_target_by_name('C_deg_A').synthesised_targets == []


def test_deg_with_interaction_as_inhibition() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_[x]_ppi+_B_[y]
                                                       A_[x]_ppi+_D_[z]
                                                       A_[x]_ppi+_E_[w]
                                                       C_deg_A; x A_[x]--B_[y]""").rxncon_system)

    for rxn in ['C_deg_A#i0', 'C_deg_A#i1']:
        target = boolean_model.reaction_target_by_name(rxn)
        assert not target.synthesised_targets
        assert target_from_str('A_[x]--0') in target.degraded_targets
        assert target_from_str('A_[x]--D_[z]') in target.degraded_targets
        assert target_from_str('A_[x]--E_[w]') in target.degraded_targets

        if not target.consumed_targets:
            assert not target.produced_targets
        elif target.consumed_targets == [target_from_str('A_[x]--D_[z]')]:
            assert target.produced_targets == [target_from_str('D_[z]--0')]
        elif target.consumed_targets == [target_from_str('A_[x]--E_[w]')]:
            assert target.produced_targets == [target_from_str('E_[w]--0')]
        else:
            raise AssertionError


def test_smooth_production_sources() -> None:
    """Testing Boolean model generation with smoothing function."""
    boolean_model = boolean_model_from_rxncon(Quick("""A_[b]_ppi+_B_[a]; ! A_[(r)]-{p}
                                                    A_[b]_ppi-_B_[a]
                                                    C_p+_A_[(r)]
                                                    D_p-_A_[(r)]""").rxncon_system,
                                              SmoothingStrategy.smooth_production_sources)

    expected = '(( A_[(r)]-{p} | A_[(r)]-{0} ) & ( A_[b]--0 | A_[b]--B_[a] )) & ' \
               '( A_[b]_ppi-_B_[a] & ( A_[b]--B_[a] | ( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 )) | ' \
               '( A_[b]--0 & ~( A_[b]_ppi+_B_[a] & A_[b]--0 & B_[a]--0 )))'

    for rule in boolean_model.update_rules:
        if rule.target == target_from_str('A_[b]--0'):
            assert rule.factor.is_equivalent_to(venn_from_str(expected, target_from_str))


def test_state_target_properties() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""C_deg_A
                                                    C_syn_A
                                                    C_p+_A_[(x)]
                                                    C_p-_A_[(x)]""").rxncon_system)

    neutral_state = target_from_str('A_[(x)]-{0}')

    assert isinstance(neutral_state, StateTarget)
    assert boolean_model.state_target_by_name('A_[(x)]-{p}').is_produced_by(
        boolean_model.reaction_target_by_name('C_p+_A_[(x)]'))
    assert boolean_model.state_target_by_name('A_[(x)]-{p}').is_consumed_by(
        boolean_model.reaction_target_by_name('C_p-_A_[(x)]'))
    assert boolean_model.state_target_by_name('A_[(x)]-{p}').is_degraded_by(
        boolean_model.reaction_target_by_name('C_deg_A'))
    assert boolean_model.reaction_target_by_name('C_syn_A').synthesises_component(neutral_state.components[0])
    assert boolean_model.state_target_by_name('A_[(x)]-{0}').is_synthesised_by(
        boolean_model.reaction_target_by_name('C_syn_A'))


def test_set_initial_condition() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""C_p+_A_[(x)]""").rxncon_system)

    neutral_state = target_from_str('A_[(x)]-{0}')
    assert boolean_model.initial_conditions.target_to_value[neutral_state] == True

    boolean_model.set_initial_condition(neutral_state, False)
    assert boolean_model.initial_conditions.target_to_value[neutral_state] == False


def test_deg_of_component_without_states() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""C_deg_A""").rxncon_system)

    target_to_factor = {rule.target: rule.factor for rule in boolean_model.update_rules}

    expected_rules = {
        'C_deg_A': '( C & A )',
        'A': '( A & ~( C_deg_A ))',
        'C': 'C'
    }
    for target_str, factor_str in expected_rules.items():
        assert target_to_factor[target_from_str(target_str)].is_equivalent_to(venn_from_str(factor_str,
                                                                                            target_from_str))


def test_boolnet_export() -> None:
    model = boolean_model_from_rxncon(Quick("""C_deg_A""").rxncon_system)
    boolnet_str, mapping, init_values = boolnet_from_boolean_model(model)
    assert 'targets, factors' in boolnet_str
    assert all(target in mapping.values() for target in
               list(model._reaction_targets.keys()) + list(model._state_targets.keys()))
    assert len(mapping) == len(init_values)
    assert all(
        init_values[key] == model.initial_conditions.target_to_value[target_from_str(value)] for key, value in
        mapping.items())


def test_homodimer_degradation() -> None:
    """The degradation of homodimers should not lead to the production of the partner, but to the complete degradation
    of the complex."""
    boolean_model = boolean_model_from_rxncon(Quick("""UC132_deg_Ste5
             Ste5_[Ste7]_ppi+_Ste7_[Ste5]
             Ste4_[Ste5]_ppi+_Ste5_[Ste4]
             Ste5_[Ste5]_ppi+_Ste5_[Ste5]
             """).rxncon_system)

    num_degs = 0

    for rule in (x for x in boolean_model.update_rules if isinstance(x.target, ReactionTarget)):
        assert isinstance(rule.target, ReactionTarget)
        if rule.target.degraded_targets:
            assert target_from_str('Ste5_[Ste5]--Ste5_[Ste5]') in rule.target.degraded_targets
            num_degs += 1
            if target_from_str('Ste5_[Ste7]--Ste7_[Ste5]') in rule.target.consumed_targets:
                assert target_from_str('Ste7_[Ste5]--0') in rule.target.produced_targets
            elif target_from_str('Ste4_[Ste5]--Ste5_[Ste4]') in rule.target.consumed_targets:
                assert target_from_str('Ste4_[Ste5]--0') in rule.target.produced_targets
            else:
                assert False

    assert num_degs == 2


def test_single_input_not_output() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_p+_B_[(a)]; ! [global]
                                                        """).rxncon_system)

    # Component expressions.
    A = 'A'
    B = '( B_[(a)]-{0} | B_[(a)]-{p} )'

    expected_rules = {
        'A_p+_B_[(a)]': '{0} & {1} & [global]'.format(A, B),
        '[global]': '[global]'
    }

    rules_found = 0

    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'A_p+_B_[(a)]':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[global]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1

    assert rules_found == 2


def test_no_input_single_output() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_p+_B_[(a)]
                                                        [global]; ! B_[(a)]-{p}""").rxncon_system)
    # Component expressions.
    A = 'A'
    B = '( B_[(a)]-{0} | B_[(a)]-{p} )'

    expected_rules = {
        'A_p+_B_[(a)]': '{0} & {1}'.format(A, B),
        '[global]': 'B_[(a)]-{p}'
    }

    rules_found = 0
    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'A_p+_B_[(a)]':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[global]':
            assert isinstance(update_rule.target, ReactionTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1

    assert rules_found == 2


def test_matching_input_output() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_p+_B_[(a)]; ! [global]
                                                        [global]; ! B_[(a)]-{p}""").rxncon_system)

    # Component expressions.
    A = 'A'
    B = '( B_[(a)]-{0} | B_[(a)]-{p} )'

    expected_rules = {
        'A_p+_B_[(a)]': '{0} & {1} & [global]'.format(A, B),
        '[global]': 'B_[(a)]-{p}'
    }

    rules_found = 0
    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'A_p+_B_[(a)]':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[global]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1

    assert rules_found == 2


def test_multiple_matching_input_one_output() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_p+_B_[(a)]; ! [global]
                                                        A_p+_B_[(a1)]; ! [global]
                                                        [global]; ! B_[(a)]-{p}""").rxncon_system)
    # Component expressions.
    A = 'A'
    B = '( B_[(a)]-{0} | B_[(a)]-{p} ) & ( B_[(a1)]-{0} | B_[(a1)]-{p} )'

    expected_rules = {
        'A_p+_B_[(a)]': '{0} & {1} & [global]'.format(A, B),
        'A_p+_B_[(a1)]': '{0} & {1} & [global]'.format(A, B),
        '[global]': 'B_[(a)]-{p}'
    }

    assert [str(update_rule.target) for update_rule in boolean_model.update_rules].count('[global]') == 1

    rules_found = 0
    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'A_p+_B_[(a)]':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == 'A_p+_B_[(a1)]':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[global]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1

    assert rules_found == 3


def test_matching_non_matching_input_one_output() -> None:
    boolean_model = boolean_model_from_rxncon(Quick("""A_p+_B_[(a)]; ! [global]
                                                        A_p+_B_[(a1)]; ! [global_diff]
                                                        [global]; ! B_[(a)]-{p}""").rxncon_system)

    # Component expressions.
    A = 'A'
    B = '( B_[(a)]-{0} | B_[(a)]-{p} ) & ( B_[(a1)]-{0} | B_[(a1)]-{p} )'

    expected_rules = {
        'A_p+_B_[(a)]': '{0} & {1} & [global]'.format(A, B),
        'A_p+_B_[(a1)]': '{0} & {1} & [global_diff]'.format(A, B),
        '[global]': 'B_[(a)]-{p}',
        '[global_diff]': '[global_diff]'
    }

    rules_found = 0

    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'A_p+_B_[(a)]':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == 'A_p+_B_[(a1)]':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[global]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[global_diff]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
    assert rules_found == 4


def test_degradation_inhibited_by_input_one_output() -> None:
    rxncon_sys = Quick("""APC_Ub+_Ndd1
                        [Ndd1UB]; ! Ndd1-{ub}
                        Decay_DEG_Ndd1mRNA; x [Ndd1UB]""").rxncon_system
    boolean_model = boolean_model_from_rxncon(rxncon_sys)

    # Component expressions.
    Ndd1mRNA = 'Ndd1mRNA'
    Decay = 'Decay'

    expected_rules = {
        'Decay_deg_Ndd1mRNA': '{0} & {1} & ~( [Ndd1UB] )'.format(Ndd1mRNA, Decay),
        '[Ndd1UB]': 'Ndd1_[(APC)]-{ub}',
    }
    rules_found = 0

    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'Decay_deg_Ndd1mRNA':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[Ndd1UB]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
    assert rules_found == 2


def test_degradation_inhibited_by_input_no_output() -> None:
    rxncon_sys = Quick("""Decay_DEG_Ndd1mRNA; x [Ndd1UB]""").rxncon_system
    boolean_model = boolean_model_from_rxncon(rxncon_sys)

    # Component expressions.
    Ndd1mRNA = 'Ndd1mRNA'
    Decay = 'Decay'

    expected_rules = {
        'Decay_deg_Ndd1mRNA': '{0} & {1} & ~( [Ndd1UB] )'.format(Ndd1mRNA, Decay),
        'Ndd1mRNA': '{0} & ~( Decay_deg_Ndd1mRNA )'.format(Ndd1mRNA),
        '[Ndd1UB]': '[Ndd1UB]',
    }
    rules_found = 0

    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'Decay_deg_Ndd1mRNA':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == 'Ndd1mRNA':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[Ndd1UB]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
    assert rules_found == 3


def test_degradation_input_required_one_output() -> None:
    rxncon_sys = Quick("""APC_Ub+_Ndd1
                        [Ndd1UB]; ! Ndd1-{ub}
                        Decay_DEG_Ndd1mRNA; ! [Ndd1UB]""").rxncon_system
    boolean_model = boolean_model_from_rxncon(rxncon_sys)

    # Component expressions.
    Ndd1mRNA = 'Ndd1mRNA'
    Decay = 'Decay'

    expected_rules = {
        'Decay_deg_Ndd1mRNA': '{0} & {1} & [Ndd1UB]'.format(Ndd1mRNA, Decay),
        'Ndd1mRNA': '{0} & ~( Decay_deg_Ndd1mRNA )'.format(Ndd1mRNA),
        '[Ndd1UB]': 'Ndd1_[(APC)]-{ub}',
    }
    rules_found = 0

    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'Decay_deg_Ndd1mRNA':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == 'Ndd1mRNA':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[Ndd1UB]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
    assert rules_found == 3


def test_degradation_input_required_no_output() -> None:
    rxncon_sys = Quick("""Decay_DEG_Ndd1mRNA; ! [Ndd1UB]""").rxncon_system
    boolean_model = boolean_model_from_rxncon(rxncon_sys)

    # Component expressions.
    Ndd1mRNA = 'Ndd1mRNA'
    Decay = 'Decay'

    expected_rules = {
        'Decay_deg_Ndd1mRNA': '{0} & {1} & [Ndd1UB]'.format(Ndd1mRNA, Decay),
        '[Ndd1UB]': '[Ndd1UB]',
    }
    rules_found = 0

    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'Decay_deg_Ndd1mRNA':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[Ndd1UB]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
    assert rules_found == 2


def test_degradation_boolean_OR_required_input_state() -> None:
    rxncon_sys = Quick("""APC_Ub+_Ndd1
                        Decay_DEG_Ndd1; ! <bool>
                        <bool>; OR [Ndd1UB]
                        <bool>; OR Ndd1-{ub}""").rxncon_system
    boolean_model = boolean_model_from_rxncon(rxncon_sys)

    # Component expressions.
    Ndd1 = '( Ndd1_[(APC)]-{0} | Ndd1_[(APC)]-{ub} )'
    Decay = 'Decay'

    target_to_factor = {rule.target: rule.factor for rule in boolean_model.update_rules}

    # C_deg_A#c0 degrades A_[(r1)]-{p}, C_deg_A#c1 degrades A_[(r2)]-{p}
    if target_to_factor[target_from_str('Decay_deg_Ndd1#c0')].is_equivalent_to(
            venn_from_str('( {0} & {1} & [Ndd1UB] )'.format(Decay, Ndd1),
                          target_from_str)):
        expected_rules = {
            'Decay_deg_Ndd1#c0': '{0} & {1} & [Ndd1UB]'.format(Ndd1, Decay),
            'Decay_deg_Ndd1#c1': '{0} & {1} & Ndd1_[(APC)]-{{ub}}'.format(Ndd1, Decay),
            'Ndd1_[(APC)]-{0}': 'Ndd1_[(APC)]-{{0}} & ~( Ndd1_[(APC)]-{{0}} & '
                                'APC_Ub+_Ndd1 ) & {0} & ~( Decay_deg_Ndd1#c0 )'.format(Ndd1),
            'Ndd1_[(APC)]-{ub}': '( Ndd1_[(APC)]-{{0}} & APC_Ub+_Ndd1 ) & ~( Decay_deg_Ndd1#c0 | Decay_deg_Ndd1#c1 ) &'
                                 ' {0} | Ndd1_[(APC)]-{{ub}} & {0} & ~( Decay_deg_Ndd1#c0 | Decay_deg_Ndd1#c1 )'
                                 .format(Ndd1),
            '[Ndd1UB]': '[Ndd1UB]',
        }

    elif target_to_factor[target_from_str('Decay_deg_Ndd1#c1')].is_equivalent_to(
            venn_from_str('( {0} & {1} & [Ndd1UB] )'.format(Decay, Ndd1),
                          target_from_str)):
        expected_rules = {
            'Decay_deg_Ndd1#c1': '{0} & {1} & [Ndd1UB]'.format(Ndd1, Decay),
            'Decay_deg_Ndd1#c0': '{0} & {1} & Ndd1_[(APC)]-{{ub}}'.format(Ndd1, Decay),
            'Ndd1_[(APC)]-{0}': 'Ndd1_[(APC)]-{{0}} & ~( Ndd1_[(APC)]-{{0}} & APC_Ub+_Ndd1 ) & '
                                '{0} & ~( Decay_deg_Ndd1#c0 )'.format(Ndd1),
            'Ndd1_[(APC)]-{ub}': '( Ndd1_[(APC)]-{{0}} & APC_Ub+_Ndd1 ) & ~( Decay_deg_Ndd1#c0 | Decay_deg_Ndd1#c1 ) & '
                                 '{0} | Ndd1_[(APC)]-{{ub}} & {0} & ~( Decay_deg_Ndd1#c0 | Decay_deg_Ndd1#c1 )'
                                 .format(Ndd1),
            '[Ndd1UB]': '[Ndd1UB]',
        }

    else:
        raise AssertionError

    rules_found = 0
    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'Decay_deg_Ndd1#c0':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == 'Decay_deg_Ndd1#c1':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == 'Ndd1_[(APC)]-{0}':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[Ndd1UB]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
    assert rules_found == 4


def test_degradation_boolean_AND_required_input_state() -> None:
    rxncon_sys = Quick("""APC_Ub+_Ndd1
                        Decay_DEG_Ndd1; ! <bool>
                        <bool>; AND [Ndd1UB]
                        <bool>; AND Ndd1-{ub}""").rxncon_system
    boolean_model = boolean_model_from_rxncon(rxncon_sys)

    # Component expressions.
    Ndd1 = '( Ndd1_[(APC)]-{0} | Ndd1_[(APC)]-{ub} )'
    Decay = 'Decay'

    expected_rules = {
        'Decay_deg_Ndd1': '{0} & {1} & Ndd1_[(APC)]-{{ub}} & [Ndd1UB]'.format(Ndd1, Decay),
        'Ndd1_[(APC)]-{0}': 'Ndd1_[(APC)]-{{0}} & ~( Ndd1_[(APC)]-{{0}} & APC_Ub+_Ndd1 )'.format(Ndd1),
        'Ndd1_[(APC)]-{ub}': '( Ndd1_[(APC)]-{{0}} & APC_Ub+_Ndd1 ) & ~( Decay_deg_Ndd1 ) & {0} | '
                             'Ndd1_[(APC)]-{{ub}} & {0} & ~( Decay_deg_Ndd1 )'.format(Ndd1),
        '[Ndd1UB]': '[Ndd1UB]',
    }
    rules_found = 0

    for update_rule in boolean_model.update_rules:
        if str(update_rule.target) == 'Decay_deg_Ndd1':
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
        elif str(update_rule.target) == '[Ndd1UB]':
            assert isinstance(update_rule.target, StateTarget)
            assert update_rule.factor.is_equivalent_to(venn_from_str(expected_rules[str(update_rule.target)],
                                                                     target_from_str))
            rules_found += 1
    assert rules_found == 2


def test_degradation_bond_to_homodimer() -> None:
    rxncon_sys = Quick("""A_[b]_ppi_B_[a]
                          B_[b]_ppi_B_[b]
                          C_deg_A; ! <x>
                          <x>; AND A@1--B@2
                          <x>; AND B@2--B@3""").rxncon_system

    boolean_model = boolean_model_from_rxncon(rxncon_sys)

    assert boolean_model.reaction_target_by_name('C_deg_A').consumed_targets == [target_from_str('A_[b]--B_[a]')]
    assert boolean_model.reaction_target_by_name('C_deg_A').produced_targets == [target_from_str('B_[a]--0')]