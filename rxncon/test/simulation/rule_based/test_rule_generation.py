from collections import namedtuple
from typing import List
import pytest
from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, rule_from_string
from rxncon.simulation.rule_based.rule_based_model import Rule
from rxncon.simulation.rule_based.rbm_from_rxncon import rules_from_reaction
from rxncon.semantics.molecule_definition_from_rxncon import mol_defs_from_rxncon_sys


RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'mol_def_strings', 'rule_strings'])

"""

LRbind: L(r,r) + R(l) <-> L(r,r!1).R(l!1) kp1, km1
LRdimer: L(r,r!1).R(l!1) + R(l) <-> L(r!2,r!1).R(l!1).R(l!2) kp2, km2

"""

def test_rule_generation(case_covalent_modifications_quant_contingencies):
    for test_case in case_covalent_modifications_quant_contingencies:
        assert is_rule_test_case_correct(test_case)

# DONE
@pytest.fixture
def case_basic_covalent_modification():
    return [
        RuleTestCase(
            'A_p+_B',
            ['A#', 'B#mod/B_[(Asite)]:u~p'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:p @ k_A_p+_B']
        ),
        # @todo fix this.
        # RuleTestCase(
        #     'A_ap_B',
        #     ['A#', 'B#mod/B_[(Asite)]:u~p'],
        #     ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:p']
        # ),
        RuleTestCase(
            'A_pt_B',
            ['A#mod/A_[(Bsite)]:u~p', 'B#mod/B_[(Asite)]:u~p'],
            ['A#mod/A_[(Bsite)]:p + B#mod/B_[(Asite)]:u -> A#mod/A_[(Bsite)]:u + B#mod/B_[(Asite)]:p @ k_A_pt_B']
        ),
        RuleTestCase(
            'A_p-_B',
            ['A#', 'B#mod/B_[(Asite)]:u~p'],
            ['A# + B#mod/B_[(Asite)]:p -> A# + B#mod/B_[(Asite)]:u @ k_A_p-_B']
        ),
        RuleTestCase(
            'A_gef_B',
            ['A#', 'B#mod/B_[(Asite)]:u~gtp'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:gtp @ k_A_gef_B']
        ),
        RuleTestCase(
            'A_gap_B',
            ['A#', 'B#mod/B_[(Asite)]:u~gtp'],
            ['A# + B#mod/B_[(Asite)]:gtp -> A# + B#mod/B_[(Asite)]:u @ k_A_gap_B']
        ),
        RuleTestCase(
            'A_ub+_B',
            ['A#', 'B#mod/B_[(Asite)]:u~ub'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:ub @ k_A_ub+_B']
        ),
        RuleTestCase(
            'A_ub-_B',
            ['A#', 'B#mod/B_[(Asite)]:u~ub'],
            ['A# + B#mod/B_[(Asite)]:ub -> A# + B#mod/B_[(Asite)]:u @ k_A_ub-_B']
        ),
        RuleTestCase(
            'A_cut_B',
            ['A#', 'B#mod/B_[(Asite)]:u~truncated'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:truncated @ k_A_cut_B']
        ),
    ]


# DONE
@pytest.fixture
def case_covalent_modifications_strict_contingencies():
    return [
        RuleTestCase(
            '''
            A_p+_B; ! A--C
            A_ppi_C''',
            ['A#ass/A_[Cassoc]:C_[Aassoc]', 'B#mod/B_[(Asite)]:u~p', 'C#ass/C_[Aassoc]:A_[Cassoc]'],
            ['A#ass/A_[Cassoc]:C_[Aassoc].C#ass/C_[Aassoc]:A_[Cassoc] + B#mod/B_[(Asite)]:u -> A#ass/A_[Cassoc]:C_[Aassoc].C#ass/C_[Aassoc]:A_[Cassoc] + B#mod/B_[(Asite)]:p @ k_A_p+_B',
             'A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#ass/A_[Cassoc]:C_[Aassoc].C#ass/C_[Aassoc]:A_[Cassoc] @ kf_A_ppi_C, kr_A_ppi_C']
        ),
        RuleTestCase(
            '''
            A_p+_B; x A--C
            A_ppi_C''',
            ['A#ass/A_[Cassoc]:C_[Aassoc]', 'B#mod/B_[(Asite)]:u~p', 'C#ass/C_[Aassoc]:A_[Cassoc]'],
            ['A#ass/A_[Cassoc]: + B#mod/B_[(Asite)]:u -> A#ass/A_[Cassoc]: + B#mod/B_[(Asite)]:p @ k_A_p+_B',
             'A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#ass/A_[Cassoc]:C_[Aassoc].C#ass/C_[Aassoc]:A_[Cassoc] @ kf_A_ppi_C, kr_A_ppi_C']
        ),
        RuleTestCase(
            '''
            A_p+_B
            B_ppi_C; ! B-{p}''',
            ['A#', 'B#ass/B_[Cassoc]:C_[Bassoc],mod/B_[(Asite)]:u~p', 'C#ass/C_[Bassoc]:B_[Cassoc]'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:p @ k_A_p+_B',
             'B#mod/B_[(Asite)]:p,ass/B_[Cassoc]: + C#ass/C_[Bassoc]: <-> B#mod/B_[(Asite)]:p,ass/B_[Cassoc]:C_[Bassoc].C#ass/C_[Bassoc]:B_[Cassoc] @ kf_B_ppi_C, kr_B_ppi_C']
        ),
        RuleTestCase(
            '''
            D_ppi_E
            A_ppi_B
            D_pt_A; ! <comp>
            <comp>; AND D--E
            <comp>; AND A--B
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(Dsite)]:u~p', 'B#ass/B_[Aassoc]:A_[Bassoc]',
             'D#ass/D_[Eassoc]:E_[Dassoc],mod/D_[(Asite)]:u~p', 'E#ass/E_[Dassoc]:D_[Eassoc]'],
            ['D#ass/D_[Eassoc]: + E#ass/E_[Dassoc]: <-> D#ass/D_[Eassoc]:E_[Dassoc].E#ass/E_[Dassoc]:D_[Eassoc] @ kr_D_ppi_E, kf_D_ppi_E',
             'A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc].B#ass/B_[Aassoc]:A_[Bassoc] @ kf_A_ppi_B, kr_A_ppi_B',
             'D#mod/D_[(Asite)]:p,ass/D_[Eassoc]:E_[Dassoc].E#ass/E_[Dassoc]:D_[Eassoc] + A#mod/A_[(Dsite)]:u,ass/A_[Bassoc]:B_[Aassoc].B#ass/B_[Aassoc]:A_[Bassoc] -> D#mod/D_[(Asite)]:u,ass/D_[Eassoc]:E_[Dassoc].E#ass/E_[Dassoc]:D_[Eassoc] + A#mod/A_[(Dsite)]:p,ass/A_[Bassoc]:B_[Aassoc].B#ass/B_[Aassoc]:A_[Bassoc] @ k_D_pt_A']
        ),
        RuleTestCase(
            '''
            D_p+_C
            A_ppi_C
            C_p+_B_[(r)]; ! C-{p}
            C_ub+_B_[(r)]; ! A--C; ! C-{P}
            ''',
            ['A#ass/A_[Cassoc]:C_[Aassoc]', 'B#mod/B_[(r)]:u~p~ub', 'C#ass/C_[Aassoc]:A_[Cassoc],mod/C_[(Dsite)]:u~p', 'D#'],
            ['D# + C#mod/C_[(Dsite)]:u -> D# + C#mod/C_[(Dsite)]:p @ k_D_p+_C',
             'A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#ass/A_[Cassoc]:C_[Aassoc].C#ass/C_[Aassoc]:A_[Cassoc] @ kf_A_ppi_C, kr_A_ppi_C',
             'C#mod/C_[(Dsite)]:p + B#mod/B_[(r)]:u -> C#mod/C_[(Dsite)]:p + B#mod/B_[(r)]:p @ k_C_p+_B_[(r)]',
             'A#ass/A_[Cassoc]:C_[Aassoc].C#mod/C_[(Dsite)]:p,ass/C_[Aassoc]:A_[Cassoc] + B#mod/B_[(r)]:u -> A#ass/A_[Cassoc]:C_[Aassoc].C#mod/C_[(Dsite)]:p,ass/C_[Aassoc]:A_[Cassoc] + B#mod/B_[(r)]:ub @ k_C_ub+_B_[(r)]']
        ),
        #
        # RuleTestCase(
        #     '''
        #     Ste5_[MEKK]_ppi_Ste11
        #     Ste5_[MEK]_ppi_Ste7
        #     Ste5_[BDSte5]_ppi_Ste5_[BDSte5]
        #     Ste11_[KD]_P+_Ste7_[(ALS359)]; ! <Ste7-5-5-11>
        #     <Ste7-5-5-11>; AND Ste5_[MEKK]--Ste11; AND Ste5_[MEK]--Ste7; AND Ste5_[BDSte5]--Ste5_[BDSte5]''',
        #     ['Ste5#ass/Ste5_[MEK]:Ste7_[Ste5assoc], ass/Ste5_[MEKK]:Ste11_[Ste5assoc], ass/Ste5_[BDSte5]:Ste5_[BDSte5]',
        #      'Ste11#ass/Ste11_[Ste5assoc]:Ste5_[MEKK]','Ste7#mod/Ste7_[(ALS359)]:u~p, ass/Ste7_[Ste5assoc]:Ste5_[MEK]'
        #      ],
        #     # Ste11#ass/Ste11_Ste5assoc]:Ste5_[MEKK], Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]
        #     #                                              ass/Ste5_[MEKK]:Ste11_[Ste5assoc]
        #     # Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5], ass/Ste5_[MEK]:Ste7_[Ste5assoc]
        #     # the first rule is new and conciders that we should have the pattern Ste5(BDSte5, MEKK).Ste5(BDSte5, MEK) and Ste5(BDSte5, MEKK, MEK).Ste5(BDSte5) to
        #     # reach the entire state space
        #     ['''Ste11#ass/Ste11_[Ste5assoc]:Ste5_[MEKK]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEKK]:Ste11_[Ste5assoc]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEK]:Ste7_[Ste5assoc]~2.Ste7#mod/Ste7_[(ALS359)]:u, ass/Ste7_[Ste5assoc]:Ste5_[MEK]~2
        #      -> Ste11#ass/Ste11_[Ste5assoc]:Ste5_[MEKK]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEKK]:Ste11_[Ste5assoc]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEK]:Ste7_[Ste5assoc]~2.Ste7#mod/Ste7_[(ALS359)]:p, ass/Ste7_[Ste5assoc]:Ste5_[MEK]~2''',
        #      '''Ste11#ass/Ste11_[Ste5assoc]:Ste5_[MEK]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEK]:Ste7_[Ste5assoc]~2, ass/Ste5_[MEKK]:Ste11_[Ste5assoc]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1.Ste7#mod/Ste7_[(ALS359)]:u, ass/Ste7_[Ste5assoc]:Ste5_[MEK]~2
        #      -> Ste11#ass/Ste11_[Ste5assoc]:Ste5_[MEK]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEK]:Ste7_[Ste5assoc]~2, ass/Ste5_[MEKK]:Ste11_[Ste5assoc]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1.Ste7#mod/Ste7_[(ALS359)]:p, ass/Ste7_[Ste5assoc]:Ste5_[MEK]~2'''
        #      ]
        # )
    ]


@pytest.fixture
def case_covalent_modifications_quant_contingencies():
    return [
        RuleTestCase(
            '''
            D_p+_C
            C_p+_B_[(r)]; k+ C-{p}''',
            ['B#mod/B_[(r)]:u~p', 'C#mod/C_[(Dsite)]:u~p', 'D#'],
            ['D# + C#mod/C_[(Dsite)]:u -> D# + C#mod/C_[(Dsite)]:p',
             'C#mod/C_[(Dsite)]:p + B#mod/B_[(r)]:u -> C#mod/C_[(Dsite)]:p + B#mod/B_[(r)]:p',
             'C#mod/C_[(Dsite)]:u + B#mod/B_[(r)]:u -> C#mod/C_[(Dsite)]:u + B#mod/B_[(r)]:p']
        ),

        RuleTestCase(
            '''
            D_p+_C
            A_ppi_C
            C_ub+_B_[(r)]; k+ A--C; ! C-{P}''',
            ['A#ass/A_[Cassoc]:C_[Aassoc]', 'B#mod/B_[(r)]:u~ub', 'C#ass/C_[Aassoc]:A_[Cassoc],mod/C_[(Dsite)]:u~p',
             'D#'],
            ['D# + C#mod/C_[(Dsite)]:u -> D# + C#mod/C_[(Dsite)]:p',
             'A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0',
             '''A#ass/A_[Cassoc]:C_[Aassoc]~0.C#mod/C_[(Dsite)]:p,ass/C_[Aassoc]:A_[Cassoc]~0 + B#mod/B_[(r)]:u
             -> A#ass/A_[Cassoc]:C_[Aassoc]~0.C#mod/C_[(Dsite)]:p,ass/C_[Aassoc]:A_[Cassoc]~0 + B#mod/B_[(r)]:ub''',
             'C#mod/C_[(Dsite)]:p,ass/C_[Aassoc]: + B#mod/B_[(r)]:u -> C#mod/C_[(Dsite)]:p,ass/C_[Aassoc]: + B#mod/B_[(r)]:ub'
             ]
        )
    ]


@pytest.fixture
def case_basic_interaction():
    return [
            RuleTestCase(
                'A_ppi_B',
                ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Aassoc]:A_[Bassoc]'],
                ['A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0']
                # we could change this to
                # A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: -> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0
                # A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0 -> A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]:
                # this would make the application of indirect contingencies a bit easier (guess) todo: discussion
                ),

           RuleTestCase(
                'A_i_B',
                ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Aassoc]:A_[Bassoc]'],
                ['A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0']
                ),

           RuleTestCase(
                'A_bind_B',
                ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Aassoc]:A_[Bassoc]'],
                ['A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0']
                ),

           RuleTestCase(
                'A_[a]_ipi_A_[b]',
                ['A#ass/A_[a]:A_[b], ass/A_[b]:A_[a]'],
                ['A#ass/A_[a]:, ass/A_[c]: <-> A#ass/A_[a]:A_[c]~0, ass/A_[c]:A_[a]~0']
                )
    ]


@pytest.fixture
def case_basic_synthesis_degradation():
    return [
        RuleTestCase(
                'A_deg_B',
                ['A#', 'B#'],
                ['A# + B# -> A#']
                # or should we set
                #['A# + B# -> A# + 0']
                ),

        RuleTestCase(
                'A_syn_B',
                ['A#', 'B#'],
                ['A# -> A# + B#']
                # or should we set
                # ['0 + A# -> A# + B#']
                ),

        # todo: trsc, trsl
        RuleTestCase(
            '''A_ppi_B
               C_p+_B
               D_ub+_B_[x]
               Y_trsl_B''',
            ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#mod/B_[(Csite)]:u~p, mod/B_[x/(Dsite)]:u~ub, ass/B_[Assoc]:A_[Bassoc]',
             'C#', 'D#', 'Y#'],
            ['A#ass/A_[Bassoc]: + B#ass/B_[Assoc]: <-> A#ass/A_[Bassoc]:B_[Assoc]~0 + B#ass/B_[Assoc]:A_[Bassoc]~0',
             'C# + B#mod/B_[(Csite)]:u -> C# + B#mod/B_[(Csite)]:p',
             'D# + B#mod/B_[x/(Dsite)]:u -> D# + B#mod/B_[x/(Dsite)]:ub',
             # to label it internally wiht mRNA is probably not a good idea because in a simulation this will acumulate
             # and the user has to include an additional equation to get rid of a molecule type we introduced.
             # TODO: discussion
             'Y# + BmRNA# -> Y# + BmRNA# + B#mod/B_[(Csite)]:u, mod/B_[x/(Dsite)]:u, ass/B_[Assoc]:']
        ),

        # todo: discuss
        RuleTestCase(
            '''
            Y_trsc_B
            ''',
            # todo: here the issue starts we transcribe B but the gene B not the protein which are two different sets.
            # The protein B can be p+ and trunc and interact, but the gene B is static (in principle).
            # What we are saying here is that we transcribe a Protein or we say later with A_p+_B that we
            # phosphorilate a gene. We have to distinguish both the protein and the gene in a better way otherwise we
            # get ambiguous

            ['B#', 'Y#', 'BmRNA#'],
            ['Y# -> Y# + BmRNA#']
        )
    ]


@pytest.fixture
def case_synthesis_degredation_contingencies():
    RuleTestCase(
                '''
                A_ppi_C
                C_ub+_B
                A_deg_B; ! A--C; ! B-{ub}
                ''',
                ['A#ass/A_[Cassoc]:C_[Aassoc]', 'B#mod/B_[(Csite)]:u~ub', 'C#ass/C_[Aassoc]:A_[Cassoc]'],
                ['A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0',
                 'C# + B#mod/B_[(Csite)]:u -> C# + B#mod/B_[(Csite)]:ub',
                 '''A#ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0 + B#mod/B_[(Csite)]:ub
                 -> A#ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0''']
                )


@pytest.fixture
def case_interaction_with_contingencies():
    return [
        RuleTestCase(
            '''
            A_ppi_B; ! A-{p}
            C_p+_A''',
            ['A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(Csite)]:u~p', 'B#ass/B_[Aassoc]:A_[Bassoc]', 'C#'],
            ['A#ass/A_[Bassoc]:, mod/A_[(Csite)]:p + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0,mod/A_[(Csite)]:p.B#ass/B_[Aassoc]:A_[Bassoc]~0',
             'C# + A#mod/A_[(Csite)]:u -> C# + A#mod/A_[(Csite)]:p']
        ),

        RuleTestCase(
            '''
            A_ppi_B_[d/s]; ! A-{p}
            D_ppi_B_[d/s]
            C_p+_A''',
            ['A#ass/A_[Bassoc]:B_[d/s],mod/A_[(Csite)]:u~p', 'B#ass/B_[d/s]:A_[Bassoc]~D_[Bassoc]', 'D#ass/D_[Bassoc]:B_[d/s]', 'C#'],
            ['A#ass/A_[Bassoc]:, mod/A_[(Csite)]:p + B#ass/B_[d/s]: <-> A#ass/A_[Bassoc]:B_[d/s]~0, mod/A_[(Csite)]:p.B#ass/B_[d/s]:A_[Bassoc]~0',
             'D#ass/D_[Bassoc]: + B#ass/B_[d/s]: <-> B#ass/B_[d/s]:D_[Bassoc]~0.D#ass/D_[Bassoc]:B_[d/s]~0',
             'C# + A#mod/A_[(Csite)]:u -> C# + A#mod/A_[(Csite)]:p']
        ),
        RuleTestCase(
            '''
            D_p+_A
            B_[a]_ipi_B_[b]
            A_ppi_B; ! <comp>
            <comp>; AND A-{P}
            <comp>; AND B_[a]--[b]
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc], mod/A_[(Dsite)]:u~p',
             'B#ass/B_[Aassoc]:A_[Bassoc], ass/B_[a]:B_[b], ass/B_[b]:B_[a]', 'D#'],
            # bound numbering has to be checked according the implementation
            ['D# + A#mod/A_[(Dsite)]:u -> D# + A#mod/A_[(Dsite)]:p',
             'B#ass/B_[a]:, ass/B_[b]: <-> B#ass/B_[a]:B_[b]~0, ass/B_[b]:B_[a]~0',
             '''A#mod/A_[(Dsite)]:p, ass/A_[Bassoc]: + B#ass/B_[Aassoc]:, ass/B_[a]:B_[b]~0, ass/B_[b]:B_[a]~0
             <-> A#mod/A_[(Dsite)]:p, ass/A_[Bassoc]:B_[Aassoc]~1.B#ass/B_[Aassoc]:A_[Bassoc]~1, ass/B_[a]:B_[b]~0, ass/B_[b]:B_[a]~0'''
             ]
        ),

        RuleTestCase(
            '''
            D_p+_A
            B_[a]_ipi_B_[b]
            A_ppi_B; ! A-{P}; x B_[a]--[b]
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc], mod/A_[(Dsite)]:u~p',
             'B#ass/B_[Aassoc]:A_[Bassoc], ass/B_[a]:B_[b], ass/B_[b]:B_[a]', 'D#'],
            ['D# + A#mod/A_[(Dsite)]:u -> D# + A#mod/A_[(Dsite)]:p',
             'B#ass/B_[a]:, ass/B_[b]: <-> B#ass/B_[a]:B_[b]~0, ass/B_[b]:B_[a]~0',
             '''A#mod/A_[(Dsite)]:p, ass/A_[Bassoc]: + B#ass/B_[Aassoc]:, ass/B_[a]:,ass/B_[b]:
             <-> A#mod/A_[(Dsite)]:p,ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[a]:,ass/B_[b]:'''
             ]
        ),

        RuleTestCase(
            '''
            C_p+_A[Gnp]
            A_ppi_B; ! A_[Gnp]-{P}; ! B_[a]--[b]
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc], mod/A_[Gnp/(Csite)]:u~p',
             'B#ass/B_[Aassoc]:A_[Bassoc], ass/B_[a]:B_[b], ass/B_[b]:B_[a]'],
            ['''A#ass/A_[Bassoc]:, mod/A_[Gnp/(Csite)]:p + B#ass/B_[Aassoc]:, ass/B_[a]:B_[b]~0, ass/B_[b]:B_[a]~0
            <-> A#ass/A_[Bassoc]:B_[Aassoc]~0, mod/A_[Gnp/(Csite)]:p.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[a]:B_[b]~1, ass/B_[b]:B_[a]~1'''
             ]
        ),

        RuleTestCase(
            '''
            A_ppi_B; ! A_[Gnp]-{P}; x B_[a]--[b]
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc], mod/A_[Gnp/(Csite)]:u~p',
             'B#ass/B_[Aassoc]:A_[Bassoc], ass/B_[a]:B_[b], ass/B_[b]:B_[a]'],
            ['''A#ass/A_[Bassoc]:, mod/A_[Gnp/(Csite)]:p + B#ass/B_[Aassoc]:, ass/B_[a]:, ass/B_[b]:
            <-> A#ass/A_[Bassoc]:B_[Aassoc]~0, mod/A_[Gnp/(Csite)]:p.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[a]:, ass/B_[b]:'''
             ]
        )
    ]


def case_disjunction():
    return [
        RuleTestCase(
            '''
            A_ppi_C
            B_ppi_D
            A_ppi_B; ! <AorC>
            <AorC>; OR A--C; OR B--D''',
            ['A#ass/A_[Cassoc]:C_[Aassoc], ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Dassoc]:D_[Bassoc], ass/B_[Assoc]:A_[Bassoc]',
             'C#ass/C_[Aassoc]:A_[Cassoc]', 'D#/ass/D_[Bassoc]:B_[Dassoc]'
             ],
            ['''A#ass/A_[Bassoc]:, ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0 + B#ass/B_[Assoc]:
            <-> A#ass/A_[Bassoc]:B_[Aassoc]~0, ass/A_[Cassoc]:C_[Aassoc]~1.B#ass/B_[Assoc]:A_[Bassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~1''',
             '''A#ass/A_[Bassoc]:, ass/A_[Cassoc]: + B#ass/B_[Assoc]:
            <-> A#ass/A_[Bassoc]:B_[Aassoc]~0, ass/A_[Cassoc]:B#ass/B_[Assoc]:A_[Bassoc]~0'''
             ]
        ),

        RuleTestCase(
            '''C_p+_A_[(x)]
               C_p+_B_[(y)]
               C_p+_B_[(z)]
               A_ppi_C
               B_ppi_D
               A_ppi_B; ! <AorC>
               <AorC>; OR A_[(x)]-{p}; OR B_[(y)]-{p}; OR A--C; OR B--D; OR B_[(z)]-{p}''',
            ['A#mod/A_[(x)]:u~p, ass/A_[Cassoc]:C_[Aassoc], ass/A_[Bassoc]:B_[Aassoc]',
             'B#mod/B_[(y)]:u~p, mod/B_[(z)]:u~p, ass/B_[Dassoc]:D_[Bassoc], ass/B_[Aassoc]:A_[Bassoc]',
             'C#ass/C_[Aassoc]:A_[Cassoc]',
             'D#ass/D_[Bassoc]:B_[Dassoc]'],
            ['C# + A#mod/A_[(x)]~u -> C# + A#mod/A_[(x)]~p',
             'C# + B#mod/B_[(y)]~u -> C# + B#mod/B_[(y)]~p',
             'C# + B#mod/B_[(z)]~u -> C# + B#mod/B_[(z)]~p',
             'A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#/ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0',
             'B#ass/B_[Dassoc]: + D#ass/D_[Bassoc]: <-> B#ass/B_[Dassoc]:D_[Bassoc]~0.D#ass/D_[Bassoc]:B_[Dassoc]~0',
             # numbering has to be checked according the implementation
              'A#mod/A_[(x)]:p, ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#mod/A_[(x)]:p, ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0',
              '''A#mod/A_[(x)]:u, ass/A_[Bassoc]:, ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0 + B#ass/B_[Aassoc]:
              <-> A#mod/A_[(x)]:u, ass/A_[Bassoc]:B_[Aassoc]~1, ass/A_[Cassoc]:C_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~1.C#ass/C_[Aassoc]:A_[Cassoc]~0''',
              '''A#mod/A_[(x)]:u, ass/A_[Bassoc]:, ass/A_[Cassoc]: + B#mod/B_[(y)]:p, ass/B_[Aassoc]:
              <-> A#mod/A_[(x)]:u, ass/A_[Bassoc]:B_[Aassoc]~0, ass/A_[Cassoc]:.B#mod/B_[(y)]:p, ass/B_[Aassoc]:A_[Bassoc]~0''',
              '''A#mod/A_[(x)]:u, ass/A_[Bassoc]:, ass/A_[Cassoc]: + B#mod/B_[(y)]:u, ass/B_[Aassoc]:, ass/B_[Dassoc]:D_[Bassoc]~0.D#ass/D_[Bassoc]:B_[Dassoc]~0
               <-> A#mod/A_[(x)]:u, ass/A_[Bassoc]:B_[Aassoc]~1, ass/A_[Cassoc]:.B#mod/B_[(y)]:u, ass/B_[Aassoc]:A_[Bassoc]~1, ass/B_[Dassoc]:D_[Bassoc]~0.D#ass/D_[Bassoc]:B_[Dassoc]~0''',
              '''A#mod/A_[(x)]:u, ass/A_[Bassoc]:, ass/A_[Cassoc]: + B#mod/B_[(y)]:u, mod/B_[(z)]:p, ass/B_[Aassoc]:, ass/B_[Dassoc]:
              <-> A#mod/A_[(x)]:u, ass/A_[Bassoc]:B_[Aassoc]~0, ass/A_[Cassoc]:.B#mod/B_[(y)]:u, mod/B_[(z)]:p, ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Dassoc]:''']
        ),

        RuleTestCase(
          '''
            A_ppi_C
            C_ppi_D
            B_ppi_E
            A_ppi_B; ! <comp1>
            <comp1>; OR <comp1C1>
            <comp1>; OR <comp2C1>
            <comp1C1>; AND A--C
            <comp1C1>; AND C--D
            <comp2C1>; AND A--C
            <comp2C1>; AND B--E''',
            ['A#ass/A_[Cassoc]:C_[Aassoc], ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Eassoc]:E_[Bassoc], ass/B_[Aassoc]:A_[Bassoc]',
             'C#ass/C_[Aassoc]:A_[Cassoc], ass/C_[Dassoc]:D_[Cassoc]', 'D#ass/D_[Cassoc]:C_[Dassoc]', 'E#ass/E_[Bassoc]:B_[Eassoc]'],
            ['A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0',
             'C#ass/C_[Dassoc]: + D#ass/D_[Cassoc]: <-> C#ass/C_[Dassoc]:D_[Cassoc]~0.D#ass/D_[Cassoc]:C_[Dassoc]~0',
             'B#ass/B_[Eassoc]: + E#ass/E_[Bassoc]: <-> B#ass/B_[Eassoc]:E_[Bassoc]~0.E#ass/E_[Bassoc]:B_[Eassoc]~0',
             # bound numbering might differ according to implementation
             #              A--C    C--D    B--E
             # [A--C, B--E]  T               T
             # [A--C, C--D]  T       T       F
             #
             '''Aass/A_[Bassoc]:, #ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0 + B#ass/B_[Aassoc]:, ass/B_[Eassoc]:E_[Bassoc]~0.E#ass/E_[Bassoc]:B_[Eassoc]~0
              <-> Aass/A_[Bassoc]:B_[Aassoc]~0, #ass/A_[Cassoc]:C_[Aassoc]~1.C#ass/C_[Aassoc]:A_[Cassoc]~1.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Eassoc]:E_[Bassoc]~2.E#ass/E_[Bassoc]:B_[Eassoc]~2''',
             '''Aass/A_[Bassoc]:, #ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0, ass/C_[Dassoc]:D_[Cassoc]~1.D#ass/D_[Cassoc]:C_[Dassoc]~1 + B#ass/B_[Aassoc]:, ass/B_[Eassoc]:
              <-> Aass/A_[Bassoc]:B_[Aassoc]~0, #ass/A_[Cassoc]:C_[Aassoc]~1.C#ass/C_[Aassoc]:A_[Cassoc]~1, ass/C_[Dassoc]:D_[Cassoc]~2.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Eassoc]:.D#ass/D_[Cassoc]:C_[Dassoc]~2'''
             ]

        ),
    ]


def case_indirect_depenendcies():
    return [
        # todo: we should think about handling ppis as two reactions instead of one lumped
        RuleTestCase(
            '''
            A_ppi_B
            B_ppi_C; ! A--B
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Aassoc]:A_[Bassoc], ass/B_[Cassoc]:C_[Bassoc]',
             'C#ass/C_[Bassoc]:B_[Cassoc]'],
            [
             # indirect dependency A_ppi_B
             # forward direction
             'A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: -> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0',
             # reverse direction k+ of B--C
             '''A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Cassoc]:C_[Bassoc]~1.C#ass/C_[Bassoc]:B_[Cassoc]~1
              -> A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]:, ass/B_[Cassoc]: + C#ass/C_[Bassoc]:''',
             '''A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Cassoc]:
              -> A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]:, ass/B_[Cassoc]:''',
             # B_ppi_C; ! A--B
             '''A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Cassoc]: + C#ass/C_[Bassoc]:
              <-> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Cassoc]:C_[Bassoc]~1.C#ass/C_[Bassoc]:B_[Cassoc]~1
             ''']
        ),

        RuleTestCase(
            '''
            Z_p+_A
            A_ppi_B; ! A-{P}
            X_p-_A''',
            ['A#ass/A_[Bassoc]:B_[Aassoc], mod/A_[(Zsite)]:u~p', 'B#ass/B_[Aassoc]:A_[Bassoc]', 'X#', 'Z#'],
            ['Z# + A#mod/A_[(Zsite)]:u -> Z# + A#mod/A_[(Zsite)]:p',
             'A#mod/A_[(Zsite)]:p, ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#mod/A_[(Zsite)]:p, ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0',
             # indirect dependency A-{P} -> k+ A--B
             'X# + A#mod/A_[(Zsite)]:p, ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0 -> X# + A#mod/A_[(Zsite)]:u, ass/A_[Bassoc]: + B#ass/B_[Aassoc]:',
             'X# + A#mod/A_[(Zsite)]:p, ass/A_[Bassoc]: -> X# + A#mod/A_[(Zsite)]:u, ass/A_[Bassoc]:'
            ]

        ),

        # todo: discuss
        RuleTestCase(
            '''
            A_ppi_B
            B_ppi_C; x A--B
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Aassoc]:A_[Bassoc], ass/B_[Cassoc]:C_[Bassoc]',
             'C#ass/C_[Bassoc]:B_[Cassoc]'],
            # A_ppi_B; k+ B--C
            ['A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]:, ass/B_[Cassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Cassoc]:',
             '''A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]:, ass/B_[Cassoc]:C_[Bassoc]~0.C#ass/C_[Bassoc]:B_[Cassoc]~0
             -> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Cassoc]: + C#ass/C_[Bassoc]:''',
             # B_ppi_C; x A--B
             'B#ass/B_[Aassoc]:, ass/B_[Cassoc]: + C#ass/C_[Bassoc]: <-> B#ass/B_[Aassoc]:, ass/B_[Cassoc]:C_[Bassoc]~0.C#ass/C_[Bassoc]:B_[Cassoc]~0'
             ]
        ),

        RuleTestCase(
            '''
            Swi4_[n]_ppi_Swi4_[m]
            Swi4_BIND_SCBG1; x Swi4_[n]--Swi4_[m]
            ''',
            ['Swi4#ass/Swi4_[n]:Swi4_[m], ass/Swi4_[m]:Swi4_[n], ass/Swi4_[SCBG1assoc]:SCBG1assoc_[Swi4assoc]'
             'SCBG1#ass/SCBG1_[Swi4assoc]:Swi4_[SCBG1assoc]'],
            # Swi4_[n]_ppi_Swi4_[m]
            ['''Swi4#ass/Swi4_[n]: + Swi4#ass/Swi4_[m]: <-> Swi4#ass/Swi4_[n]:Swi4_[m]~0.Swi4#assSwi4_[m]:Swi4_[n]~0''',
             '''SCBG1#ass/SCBG1_[Swi4assoc]:Swi4_[SCBG1assoc]~0.Swi4#ass/Swi4_[SCBG1assoc]:SCBG1assoc_[Swi4assoc]~0
            + SCBG1#ass/SCBG1_[Swi4assoc]:Swi4_[SCBG1assoc]~0.Swi4#ass/Swi4_[SCBG1assoc]:SCBG1assoc_[Swi4assoc]~0
            -> Swi4#ass/Swi4_[n]:Swi4_[m]~0.Swi4#ass/Swi4_[m]:Swi4_[n]~0 + SCBG1#ass/SCBG1_[Swi4assoc]: + SCBG1#ass/SCBG1_[Swi4assoc]:''',
             # Swi4_BIND_SCBG1; x Swi4_[n]--Swi4_[m]
             '''SCBG1#ass/SCBG1_[Swi4assoc]: + Swi4#ass/Swi4_[SCBG1assoc]:, ass/Swi4_[n]:, ass/Swi4_[m]:
             <-> SCBG1#ass/SCBG1_[Swi4assoc]:Swi4_[SCBG1assoc]~0.Swi4#ass/Swi4_[SCBG1assoc]:SCBG1_[Swi4assoc]~0, ass/Swi4_[n]:, ass/Swi4_[m]:''']

        ),

        RuleTestCase(
            """
            Swi6_[c]_ppi_Swi4_[c]
            Swi4_[n]_ipi_Swi4_[c]; x Swi6_[c]--Swi4_[c]
            """,
            ['Swi4#ass/Swi4_[n]:Swi4_[c], ass/Swi4_[c]:Swi4_[n]~Swi6_[c]',
             'Swi6#ass/Swi6_[c]:Swi4_[c]'],
            # Swi6_[c]_ppi_Swi4_[c]
            ['''Swi6#ass/Swi6_[c]: + Swi4#ass/Swi4_[c]:, ass/Swi4_[n]: <-> Swi4#ass/Swi4_[c]:Swi6_[c]~0, ass/Swi4_[n]:.Swi6#ass/Swi6_[c]:Swi4_[c]~0''',
             '''Swi6#ass/Swi6_[c]: + Swi4#ass/Swi4_[c]:Swi4_[n]~0, ass/Swi4_[n]:Swi4_[c]~0
            <-> Swi4#ass/Swi4_[c]:Swi6_[c]~0, ass/Swi4_[n]:.Swi6#ass/Swi6_[c]:Swi4_[c]~0''',
             # Swi4_[n]_ipi_Swi4_[c]; x Swi6_[c]--Swi4_[c]
             'Swi4#ass/Swi4_[c]:, ass/Swi4_[n]: <-> Swi4#ass/Swi4_[c]:Swi4_[n]~0, ass/Swi4_[n]:Swi4_[c]~0']
        )

# ### simple conflict chain
#                 '''X_p-_A
#                 A_ppi_B; ! A_[X]-{P}
#                 B_ppi_C; ! A--B
#                 C_ppi_D; ! B--C'''
#     'X + A(X~P,AssocB!3).B(AssocA!3,AssocC!2).C(AssocB!2,AssocD!1).D(AssocC!1) -> X + A(X~U,AssocB) + B(AssocA,AssocC) + C(AssocB,AssocD) + D(AssocC)',
#     'X + A(X~P,AssocB) -> X + A(X~U,AssocB)',
#     'X + A(X~P,AssocB!1).B(AssocA!1,AssocC) -> X + A(X~U,AssocB) + B(AssocA,AssocC)',
#     'X + A(X~P,AssocB!2).B(AssocA!2,AssocC!1).C(AssocB!1,AssocD) -> X + A(X~U,AssocB) + B(AssocA,AssocC) + C(AssocB,AssocD)'],

# ### conflict chain with two alternative paths
#             '''X_p-_A
#                A_ppi_B; ! A_[X]-{P}
#                A_ppi_F; ! A--B
#                B_ppi_C; ! A--B
#                C_ppi_D; ! B--C'''
#     'X + A(X~P,AssocB,AssocF) -> X + A(X~U,AssocB,AssocF)',
#     'X + A(X~P,AssocB!1,AssocF).B(AssocA!1,AssocC) -> X + A(X~U,AssocB,AssocF) + B(AssocA,AssocC)',
#     'X + A(X~P,AssocB!2,AssocF!1).B(AssocA!2,AssocC).F(AssocA!1) -> X + A(X~U,AssocB,AssocF) + F(AssocA) + B(AssocA,AssocC)',
#     'X + A(X~P,AssocB!2,AssocF).B(AssocA!2,AssocC!1).C(AssocB!1,AssocD) -> X + A(X~U,AssocB,AssocF) + B(AssocA,AssocC) + C(AssocB,AssocD)',
#     'X + A(X~P,AssocB!3,AssocF).B(AssocA!3,AssocC!2).C(AssocB!2,AssocD!1).D(AssocC!1) -> X + A(X~U,AssocB,AssocF) + B(AssocA,AssocC) + C(AssocB,AssocD) + D(AssocC)',
#     'X + A(X~P,AssocB!3,AssocF!1).B(AssocA!3,AssocC!2).C(AssocB!2,AssocD).F(AssocA!1) -> X + A(X~U,AssocB,AssocF) + F(AssocA) + B(AssocA,AssocC) + C(AssocB,AssocD)',
#     'X + A(X~P,AssocB!4,AssocF!1).B(AssocA!4,AssocC!3).C(AssocB!3,AssocD!2).D(AssocC!2).F(AssocA!1) -> X + A(X~U,AssocB,AssocF) + F(AssocA) + B(AssocA,AssocC) + C(AssocB,AssocD) + D(AssocC)'],
#
    ]


def is_rule_test_case_correct(test_case) -> bool:
    rxncon = Quick(test_case.quick_string).rxncon_system

    actual_mol_defs = set(mol_defs_from_rxncon_sys(rxncon).values())
    actual_rules = set()

    for reaction in rxncon.reactions:
        actual_rules = actual_rules.union(rules_from_reaction(rxncon, reaction))

    expected_mol_defs = {mol_def_from_string(x) for x in test_case.mol_def_strings}
    expected_rules    = {rule_from_string(expected_mol_defs, x) for x in test_case.rule_strings}

    correct_mol_defs = actual_mol_defs == expected_mol_defs
    correct_rules = actual_rules == expected_rules

    if not correct_mol_defs:
        print('Expected molecule definitions:')
        print(expected_mol_defs)
        print()
        print('Actual molecule definitions:')
        print(actual_mol_defs)

    if not correct_rules:
        print('Expected rules:')
        for rule in expected_rules:
            print(rule)
            print()
        print()
        print('Actual rules:')
        for rule in actual_rules:
            print(rule)
            print()

    return correct_mol_defs and correct_rules

