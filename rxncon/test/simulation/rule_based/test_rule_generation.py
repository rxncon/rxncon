from collections import namedtuple
from typing import List
import pytest
from rxncon.input.quick.quick import Quick
from rxncon.simulation.rule_based.molecule_from_string import mol_def_from_string, rule_from_string
from rxncon.simulation.rule_based.rbm_from_rxncon import RuleBasedModelSupervisor
from rxncon.simulation.rule_based.rule_based_model import Rule


RuleTestCase = namedtuple('RuleTestCase', ['quick_string', 'mol_def_strings', 'rule_strings'])

"""

LRbind: L(r,r) + R(l) <-> L(r,r!1).R(l!1) kp1, km1
LRdimer: L(r,r!1).R(l!1) + R(l) <-> L(r!2,r!1).R(l!1).R(l!2) kp2, km2

"""

def test_rule_generation(test_cases):
    for test_case in test_cases:
        assert is_rule_test_case_correct(test_case)


def test_basic_covalent_mod_cases():
    return [
        RuleTestCase(
            'A_p+_B',
            ['A#', 'B#mod/B_[(Asite)]:u~p'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:p']
        ),

        RuleTestCase(
            'A_ap_B',
            ['A#', 'B#mod/B_[(Asite)]:u~p'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:p']
        ),

        RuleTestCase(
            'A_pt_B',
            ['A#', 'B#mod/B_[(Asite)]:u~p'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:p']
        ),

        RuleTestCase(
            'A_p-_B',
            ['A#', 'B#mod/B_[(Asite)]:u~p'],
            ['A# + B#mod/B_[(Asite)]:p -> A# + B#mod/B_[(Asite)]:u']
        ),

        RuleTestCase(
            'A_gef_B',
            ['A#', 'B#mod/B_[(Asite)]:u~gtp'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:gtp']
        ),

        RuleTestCase(
            'A_gap_B',
            ['A#', 'B#mod/B_[(Asite)]:u~gtp'],
            ['A# + B#mod/B_[(Asite)]:gtp -> A# + B#mod/B_[(Asite)]:u']
        ),

        RuleTestCase(
            'A_ub+_B',
            ['A#', 'B#mod/B_[(Asite)]:u~ub'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:ub']
        ),

        RuleTestCase(
            'A_ub-_B',
            ['A#', 'B#mod/B_[(Asite)]:u~ub'],
            ['A# + B#mod/B_[(Asite)]:ub -> A# + B#mod/B_[(Asite)]:u']
        ),

        RuleTestCase(
            'A_cut_B',
            ['A#', 'B#mod/B_[(Asite)]:u~truncated'],
            ['A# + B#mod/B_[(Asite)]:u -> A# + B#mod/B_[(Asite)]:truncated']
        ),
    ]

def test_covalent_modifications_contingencies():
    return [
        RuleTestCase(
            '''
            D_ppi_E
            A_ppi_B
            D_pt_A; ! <comp>
            <comp>; AND D--E
            <comp>; AND A--B
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(Dsite)]:u~p', 'B#ass/B_[Aassoc]:A_[Bassoc]',
             'D#ass/D_[Eassoc]:E_[Dassoc],mod/D_[(Asite)]~u~p, E#ass/E_[Dassoc]:D_[Eassoc]'],
            ['D#ass/D_[Eassoc]: + E#ass/E_[Dassoc]: <-> D#ass/D_[Eassoc]:E_[Dassoc]~0.E#ass/E_[Dassoc]:D_[Eassoc]~0',
             'A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0',
             '''D#mod/D_[(Asite)]~p,ass/D_[Eassoc]:E_[Dassoc]~0.E#ass/E_[Dassoc]:D_[Eassoc]~0 + A#mod/A_[(Dsite)]:u,ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0
             -> D#mod/D_[(Asite)]~u,ass/D_[Eassoc]:E_[Dassoc]~0.E#ass/E_[Dassoc]:D_[Eassoc]~0 + A#mod/A_[(Dsite)]:p,ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0'''
             ]
        ),

        RuleTestCase(
            '''
            D_p+_C
            A_ppi_C
            C_p+_B_[(r)]; ! C-{p}
            C_ub+_B_[(r)]; ! A--C; ! C-{P}
            ''',
            ['A#ass/A_[Cassoc]:C_[Aassoc]', 'B#mod/B_[(r)]:u~p~ub', 'C#ass/C_[Aassoc]:A_[Cassoc],mod/C_[(Dsite)]:u~p', 'D#'],
            ['D# + C#mod/C_[(Dsite)]:u -> D# + C#mod/C_[(Dsite)]:p',
             'A#ass/A_[Cassoc]: + C#ass/C_[Aassoc]: <-> A#ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0',
             'C#mod/C_[(Dsite)]:p + B#mod/B_[(r)]:u -> C#mod/C_[(Dsite)]:p + B#mod/B_[(r)]:p',
             '''A#ass/A_[Cassoc]:C_[Aassoc]~0.C#mod/C_[(Dsite)]:p,ass/C_[Aassoc]:A_[Cassoc]~0 + B#mod/B_[(r)]:u
             -> A#ass/A_[Cassoc]:C_[Aassoc]~0.C#mod/C_[(Dsite)]:p,ass/C_[Aassoc]:A_[Cassoc]~0 + B#mod/B_[(r)]:p''']
        ),
        RuleTestCase(
            '''
            Ste5_[MEKK]_ppi_Ste11
            Ste5_[MEK]_ppi_Ste7
            Ste5_[BDSte5]_ppi_Ste5_[BDSte5]
            Ste11_[KD]_P+_Ste7_[(ALS359)]; ! <Ste7-5-5-11>
            <Ste7-5-5-11>; AND Ste5_[MEKK]--Ste11; AND Ste5_[MEK]--Ste7; AND Ste5_[BDSte5]--Ste5_[BDSte5]''',
            ['Ste5#ass/Ste5_[MEK]:Ste7_[Ste5assoc], ass/Ste5_[MEK]ass/Ste5_[MEKK]:Ste11_[Ste5assoc], ass/Ste5_[BDSte5]:Ste5_[BDSte5]',
             'Ste11#ass/Ste11_[Ste5assoc]:Ste5_[MEKK]','Ste7#mod/Ste7_[(ALS359)]:u~p, ass/Ste7_[Ste5assoc]:Ste5_[MEK]'],
            # Ste11#ass/Ste11_Ste5assoc]:Ste5_[MEKK], Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]
            #                                              ass/Ste5_[MEKK]:Ste11_[Ste5assoc]
            # Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5], ass/Ste5_[MEK]:Ste7_[Ste5assoc]
            # the first rule is new and conciders that we should have the pattern Ste5(BDSte5, MEKK).Ste5(BDSte5, MEK) and Ste5(BDSte5, MEKK, MEK).Ste5(BDSte5) to
            # reach the entire state space
            ['''Ste11#ass/Ste11_Ste5assoc]:Ste5_[MEKK]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEKK]:Ste11_[Ste5assoc]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEK]:Ste7_[Ste5assoc]~2.Ste7#mod/Ste7_[(ALS359)]:u, ass/Ste7_[Ste5assoc]:Ste5_[MEK]~2
             -> Ste11#ass/Ste11_Ste5assoc]:Ste5_[MEKK]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEKK]:Ste11_[Ste5assoc]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEK]:Ste7_[Ste5assoc]~2.Ste7#mod/Ste7_[(ALS359)]:p, ass/Ste7_[Ste5assoc]:Ste5_[MEK]~2''',
             '''Ste11#ass/Ste11_[Ste5assoc]:Ste5_[MEK]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEK]:Ste7_[Ste5assoc]~2, ass/Ste5_[MEKK]:Ste11_[Ste5assoc]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1.Ste7#mod/Ste7_[(ALS359)]:u, ass/Ste7_[Ste5assoc]:Ste5_[MEK]~2
             -> Ste11#ass/Ste11_[Ste5assoc]:Ste5_[MEK]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1, ass/Ste5_[MEK]:Ste7_[Ste5assoc]~2, ass/Ste5_[MEKK]:Ste11_[Ste5assoc]~0.Ste5#ass/Ste5_[BDSte5]:Ste5_[BDSte5]~1.Ste7#mod/Ste7_[(ALS359)]:p, ass/Ste7_[Ste5assoc]:Ste5_[MEK]~2''']
        )


    ]

def test_basic_interaction_cases():
       return [
               RuleTestCase(
                'A_ppi_B',
                ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Aassoc]:A_[Bassoc]'],
                ['A#ass/A_[Bassoc]: + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0']
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

def test_basic_synthesis_degradation_cases():
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
        # RuleTestCase(
        #     '''A_ppi_B
        #        C_p+_B
        #        D_ub+_B_[x]
        #        Y_trsl_B''',
        #     ['A#ass/A_[Bassoc]:B_[Aassoc]', 'B#mod/B_[Csite]:u~p, mod/B_[x]:u~ub, ass/B_[Assoc]:A_[Bassoc]', 'C#', 'D#', 'Y#'],
        #     ['A#ass/A_[Bassoc]: + B#ass/B_[Assoc]: <-> A#ass/A_[Bassoc]:B_[Assoc]~0 + B#ass/B_[Assoc]:A_[Bassoc]~0',
        #      'C# + B#mod/B_[Csite]:u -> C# + B#mod/B_[Csite]:p',
        #      'D# + B#mod/B_[x]:u -> D# + B#mod/B_[x]:ub',
        #      'Y# + BmRNA# -> Y# + BmRNA# + B#mod/B_[Csite]:u, mod/B_[x]:u, ass/B_[Assoc]:']
        # )

    ]

def test_synthesis_degredation_contingencies():
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
#@pytest.fixture
def test_interaction_with_contingencies():
    return [
        RuleTestCase(
            '''
            A_ppi_B; ! A-{p}
            C_p+_A''',
            ['A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(Csite)]:u~p', 'B#ass/B_[Aassoc]:A_[Bassoc]', 'C#'],
            ['A#ass/A_[Bassoc]:,mod/A_[(Csite)]:p + B#ass/B_[Aassoc]: <-> A#ass/A_[Bassoc]:B_[Aassoc]~0,mod/A_[(Csite)]:p.B#ass/B_[Aassoc]:A_[Bassoc]~0',
             'C# + A#mod/A_[(Csite)]:u -> C# + A#mod/A_[(Csite)]:p']
        ),

        RuleTestCase(
            '''
            A_ppi_B_[d/s]; ! A-{p}
            D_ppi_B_[d/s]
            C_p+_A''',
            ['A#ass/A_[Bassoc]:B_[d/s],mod/A_[(Csite)]:u~p', 'B#ass/B_[d/s]:A_[Bassoc]~D_[Bassoc]', 'D#ass/D_[Bassoc]:B_[d/s]', 'C#'],
            ['A#ass/A_[Bassoc]:,mod/A_[(Csite)]:p + B#ass/B_[d/s]: <-> A#ass/A_[Bassoc]:B_[d/s]~0,mod/A_[(Csite)]:p.B#ass/B_[d/s]:A_[Bassoc]~0',
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
            ['A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(Dsite)]:u~p',
             'B#ass/B_[Aassoc]:A_[Bassoc],ass/B_[a]:B_[b],ass/B_[b]:B_[a]', 'D#'],
            # bound numbering has to be checked according the implementation
            ['D# + A#mod/A_[(Dsite)]:u -> D# + A#mod/A_[(Dsite)]:p',
             'B#ass/B_[a]:,ass/B_[b]: <-> B#ass/B_[a]:B_[b]~0,ass/B_[b]:B_[a]~0',
             '''A#mod/A_[(Dsite)]:p,ass/A_[Bassoc]: + B#ass/B_[Aassoc]:,ass/B_[a]:B_[b]~0,ass/B_[b]:B_[a]~0
             <-> A#mod/A_[(Dsite)]:p,ass/A_[Bassoc]:B_[Aassoc]~1.B#ass/B_[Aassoc]:A_[Bassoc]~1,ass/B_[a]:B_[b]~0,ass/B_[b]:B_[a]~0'''
             ]
        ),

        RuleTestCase(
            '''
            D_p+_A
            B_[a]_ipi_B_[b]
            A_ppi_B; ! A-{P}; x B_[a]--[b]
            ''',
            ['A#ass/A_[Bassoc]:B_[Aassoc],mod/A_[(Dsite)]:u~p',
             'B#ass/B_[Aassoc]:A_[Bassoc],ass/B_[a]:B_[b],ass/B_[b]:B_[a]', 'D#'],
            ['D# + A#mod/A_[(Dsite)]:u -> D# + A#mod/A_[(Dsite)]:p',
             'B#ass/B_[a]:,ass/B_[b]: <-> B#ass/B_[a]:B_[b]~0,ass/B_[b]:B_[a]~0',
             '''A#mod/A_[(Dsite)]:p,ass/A_[Bassoc]: + B#ass/B_[Aassoc]:,ass/B_[a]:,ass/B_[b]:
             <-> A#mod/A_[(Dsite)]:p,ass/A_[Bassoc]:B_[Aassoc]~0.B#ass/B_[Aassoc]:A_[Bassoc]~0,ass/B_[a]:,ass/B_[b]:'''
             ]
        ),

        RuleTestCase(
            '''
            A_ppi_C
            B_ppi_D
            A_ppi_B; ! <AorC>
            <AorC>; OR A--C; OR B--D''',
            ['A#ass/A_[Cassoc]:C_[Aassoc], ass/A_[Bassoc]:B_[Aassoc]', 'B#ass/B_[Dassoc]:D_[Bassoc], ass/B_[Assoc]:A_[Bassoc]',
             'C#ass/C_[Aassoc]:A_[Cassoc]', 'D#/ass/D_[Bassoc]:B_[Dassoc]'],
            ['''A#ass/A_[Bassoc]:, ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0 + B#ass/B_[Assoc]:
            <-> A#ass/A_[Bassoc]:B_[Aassoc]~0, ass/A_[Cassoc]:C_[Aassoc]~1.B#ass/B_[Assoc]:A_[Bassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~1''']

        #'A(AssocB,AssocC!1).C(AssocA!1) + B(AssocA) <-> A(AssocB!2,AssocC!1).B(AssocA!2).C(AssocA!1)',
        #'A(AssocB,AssocC) + B(AssocA,AssocD!1).D(AssocB!1) <-> A(AssocB!2,AssocC).B(AssocA!2,AssocD!1).D(AssocB!1)'],

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
             '''Aass/A_[Bassoc]:, #ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0 + B#ass/B_[Aassoc]:, ass/B_[Eassoc]:E_[Bassoc]~0.E#ass/E_[Bassoc]:B_[Eassoc]~0
              <-> Aass/A_[Bassoc]:B_[Aassoc]~0, #ass/A_[Cassoc]:C_[Aassoc]~1.C#ass/C_[Aassoc]:A_[Cassoc]~1.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Eassoc]:E_[Bassoc]~2.E#ass/E_[Bassoc]:B_[Eassoc]~2''',
             '''Aass/A_[Bassoc]:, #ass/A_[Cassoc]:C_[Aassoc]~0.C#ass/C_[Aassoc]:A_[Cassoc]~0, ass/C_[Dassoc]:D_[Cassoc]~1.D#ass/D_[Cassoc]:C_[Dassoc]~1 + B#ass/B_[Aassoc]:, ass/B_[Eassoc]:
              <-> Aass/A_[Bassoc]:B_[Aassoc]~0, #ass/A_[Cassoc]:C_[Aassoc]~1.C#ass/C_[Aassoc]:A_[Cassoc]~1, ass/C_[Dassoc]:D_[Cassoc]~2.B#ass/B_[Aassoc]:A_[Bassoc]~0, ass/B_[Eassoc]:.D#ass/D_[Cassoc]:C_[Dassoc]~2'''
             ]

        )


    ]


def is_rule_test_case_correct(test_case) -> bool:
    rxncon = Quick(test_case.quick_string).rxncon_system
    actual_mol_defs = set(RuleBasedModelSupervisor(rxncon).mol_defs.values())
    actual_rules    = RuleBasedModelSupervisor(rxncon).rules

    expected_mol_defs = {mol_def_from_string(x) for x in test_case.mol_def_strings}
    expected_rules    = [rule_from_string(expected_mol_defs, x) for x in test_case.rule_strings]

    return actual_mol_defs == expected_mol_defs and are_rule_lists_equivalent(actual_rules, expected_rules)


def are_rule_lists_equivalent(first_list: List[Rule], second_list: List[Rule]) -> bool:
    while first_list:
        first_rule = first_list.pop()
        for second_rule in second_list:
            if are_rules_equivalent(first_rule, second_rule):
                second_list.remove(second_rule)

    return len(second_list) == 0


def are_rules_equivalent(first_rule: Rule, second_rule: Rule) -> bool:
    # @todo we disregard the rates in this equivalence.
    return set(first_rule.left_hand_side) == set(second_rule.left_hand_side) and \
        set(first_rule.right_hand_side) == set(second_rule.right_hand_side) and \
        first_rule.arrow_type == second_rule.arrow_type