#from rxncon.core.reaction import *
from collections import namedtuple
import pytest

import rxncon.core.reaction as rxn
import rxncon.core.state as sta
import rxncon.core.specification as spec


DefinitionTestCase = namedtuple('DefinitionTestCase', ['definition', 'expected_name', 'expected_reactant_defs_pre',
                                                       'expected_reactant_defs_post', 'expected_reactants_defs',
                                                       'expected_representation_def', 'expected_variable_x',
                                                       'expected_variable_y'])

@pytest.fixture
def reaction_definitions():
    return {'phosphorylation': rxn.ReactionDef('phosphorylation',
                                                      '$x_p+_$y',
                                               {
                                                          '$x': (spec.ProteinSpec, spec.SpecificationResolution.component),
                                                          '$y': (spec.ProteinSpec, spec.SpecificationResolution.residue)
                                                      },
                                                      '$x# + $y#$y-{0} -> $x# + $y#$y-{p}'
                                               ),

            'phosphotransfer': rxn.ReactionDef('phosphotransfer',
                                                      '$x_pt_$y',
                                               {
                                                        '$x': (spec.ProteinSpec, spec.SpecificationResolution.residue),
                                                        '$y': (spec.ProteinSpec, spec.SpecificationResolution.residue)
                                                      },
                                                      '$x#$x-{p} + $y#$y-{0} -> $x#$x-{0} + $y#$y-{p}'
                                               ),
            'protein-protein-interaction': rxn.ReactionDef('protein-protein-interaction',
                                                                  '$x_ppi_$y',
                                                           {
                                                                    '$x': (spec.ProteinSpec, spec.SpecificationResolution.domain),
                                                                   '$y': (spec.ProteinSpec, spec.SpecificationResolution.domain)
                                                                  },
                                                                  '$x#$x--0 + $y#$y--0 <-> $x#$x--$y + $y#$x--$y'
                                                           ),

            'intra-protein-interaction': rxn.ReactionDef('intra-protein-interaction',
                                                            '$x_ipi_$y',
                                                         {
                                                                '$x': (spec.ProteinSpec, spec.SpecificationResolution.domain),
                                                                '$y': (spec.ProteinSpec, spec.SpecificationResolution.domain)
                                                            },
                                                            '$x#$x--0,$y--0 -> $x#$x--$y.domain'
                                                         ),

            'transcription': rxn.ReactionDef('transcription',
                                                   '$x_trsc_$y',
                                             {
                                                        '$x': (spec.ProteinSpec, spec.SpecificationResolution.component),
                                                        '$y': (spec.DnaSpec, spec.SpecificationResolution.component)
                                                    },
                                                   '$x# + $y# -> $x# + $y# + $y.mRNA#'
                                             ),

            'translation': rxn.ReactionDef('translation',
                                                  '$x_trsl_$y',
                                           {
                                                    '$x': (spec.ProteinSpec, spec.SpecificationResolution.component),
                                                    '$y': (spec.RnaSpec, spec.SpecificationResolution.component)
                                                  },
                                                  '$x# + $y# -> $x# + $y# + $y.protein#'
                                           )
            }

@pytest.fixture
def definition_test_case(reaction_definitions):
    return [
        DefinitionTestCase(reaction_definitions['phosphorylation'],
                            'phosphorylation',
                           ['$x#', '$y#$y-{0}'],
                           ['$x#', '$y#$y-{p}'],
                            '$x# + $y#$y-{0} -> $x# + $y#$y-{p}',
                            '$x_p+_$y',
                           (spec.ProteinSpec, spec.SpecificationResolution.component),
                           (spec.ProteinSpec, spec.SpecificationResolution.residue)),

        DefinitionTestCase(reaction_definitions['phosphotransfer'],
                            'phosphotransfer',
                           ['$x#$x-{p}', '$y#$y-{0}'],
                           ['$x#$x-{0}', '$y#$y-{p}'],
                            '$x#$x-{p} + $y#$y-{0} -> $x#$x-{0} + $y#$y-{p}',
                            '$x_pt_$y',
                           (spec.ProteinSpec, spec.SpecificationResolution.residue),
                           (spec.ProteinSpec, spec.SpecificationResolution.residue)),

    DefinitionTestCase(reaction_definitions['protein-protein-interaction'],
                        'protein-protein-interaction',
                       ['$x#$x--0', '$y#$y--0'],
                       ['$x#$x--$y', '$y#$x--$y'],
                        '$x#$x--0 + $y#$y--0 <-> $x#$x--$y + $y#$x--$y',
                        '$x_ppi_$y',
                       (spec.ProteinSpec, spec.SpecificationResolution.domain),
                       (spec.ProteinSpec, spec.SpecificationResolution.domain)),

    DefinitionTestCase(reaction_definitions['transcription'],
                        'transcription',
                       ['$x#', '$y#'],
                       ['$x#', '$y#', '$y.mRNA#'],
                        '$x# + $y# -> $x# + $y# + $y.mRNA#',
                        '$x_trsc_$y',
                       (spec.ProteinSpec, spec.SpecificationResolution.component),
                       (spec.DnaSpec, spec.SpecificationResolution.component)
                       ),
    DefinitionTestCase(reaction_definitions['translation'],
                        'translation',
                       ['$x#', '$y#'],
                       ['$x#', '$y#', '$y.protein#'],
                        '$x# + $y# -> $x# + $y# + $y.protein#',
                        '$x_trsl_$y',
                       (spec.ProteinSpec, spec.SpecificationResolution.component),
                       (spec.RnaSpec, spec.SpecificationResolution.component)
                       )
    ]


def test_definition(definition_test_case):
    for the_case in definition_test_case:
        _is_definition_correct(the_case)


def _is_definition_correct(the_case):
    assert the_case.definition.name == the_case.expected_name
    assert the_case.definition.reactant_defs_post == the_case.expected_reactant_defs_post
    assert the_case.definition.reactant_defs_pre == the_case.expected_reactant_defs_pre
    assert the_case.definition.reactants_defs == the_case.expected_reactants_defs
    assert the_case.definition.representation_def == the_case.expected_representation_def
    assert the_case.definition.variables_def['$x'] == the_case.expected_variable_x
    assert the_case.definition.variables_def['$y'] == the_case.expected_variable_y

#
ReactionTestCase = namedtuple('ReactionTestCase' , ['reaction_str', 'expected_definition', 'expected_reaction',
                                                    'expected_reactants_pre', 'expected_reactants_post'])


@pytest.fixture
def reaction_test_case(reaction_definitions):
    return [
        ReactionTestCase('A@0_p+_B@1_[(y)]',
                         reaction_definitions['phosphorylation'],
                         rxn.Reaction(reaction_definitions['phosphorylation'],
                                      {'$x': rxn.specification_from_string('A@0'),
                                       '$y': rxn.specification_from_string('B@1_[(y)]')
                                       }),

                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       []),
                          rxn.Reactant(spec.ProteinSpec('B', 1, spec.Domain(None, None, None)),
                                       [sta.state_from_string('B@1_[(y)]-{0}')])],

                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       []),
                          rxn.Reactant(spec.ProteinSpec('B', 1, spec.Domain(None, None, None)),
                                       [sta.state_from_string('B@1_[(y)]-{P}')])],

                         ),

        ReactionTestCase('A@0_[(x)]_pt_B@1_[(y)]',
                         reaction_definitions['phosphotransfer'],
                         rxn.Reaction(reaction_definitions['phosphotransfer'],
                                      {'$x': rxn.specification_from_string('A@0_[(x)]'),
                                       '$y': rxn.specification_from_string('B@1_[(y)]')
                                       }),
                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       [sta.state_from_string('A@0_[(x)]-{P}')]),
                          rxn.Reactant(spec.ProteinSpec('B', 1, spec.Domain(None, None, None)),
                                       [sta.state_from_string('B@1_[(y)]-{0}')])],

                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       [sta.state_from_string('A@0_[(x)]-{0}')]),
                          rxn.Reactant(spec.ProteinSpec('B', 1, spec.Domain(None, None, None)),
                                       [sta.state_from_string('B@1_[(y)]-{P}')])],
                         ),

        ReactionTestCase('A@0_[b]_ppi_B@1_[a]',
                         reaction_definitions['protein-protein-interaction'],
                         rxn.Reaction(reaction_definitions['protein-protein-interaction'],
                                      {'$x': rxn.specification_from_string('A@0_[b]'),
                                       '$y': rxn.specification_from_string('B@1_[a]')
                                       }),
                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       [sta.state_from_string('A@0_[b]--0')]),
                          rxn.Reactant(spec.ProteinSpec('B', 1, spec.Domain(None, None, None)),
                                       [sta.state_from_string('B@1_[a]--0')])],

                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       [sta.state_from_string('A@0_[b]--B@1_[a]')]),
                          rxn.Reactant(spec.ProteinSpec('B', 1, spec.Domain(None, None, None)),
                                       [sta.state_from_string('A@0_[b]--B@1_[a]')])],

                         ),
        ReactionTestCase('A@0_[n]_ipi_A@0_[m]',
                         reaction_definitions['intra-protein-interaction'],
                         rxn.Reaction(reaction_definitions['intra-protein-interaction'],
                                      {'$x': rxn.specification_from_string('A@0_[n]'),
                                       '$y': rxn.specification_from_string('A@0_[m]')
                                       }),
                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       [sta.state_from_string('A@0_[n]--0'), sta.state_from_string('A@0_[m]--0')])],

                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       [sta.state_from_string('A@0_[n]--[m]')])],

                         ),

        ReactionTestCase('A@0_trsc_Bgene@1',
                         reaction_definitions['transcription'],
                         rxn.Reaction(reaction_definitions['transcription'],
                                      {'$x': rxn.specification_from_string('A@0'),
                                       '$y': rxn.specification_from_string('Bgene@1')
                                       }),
                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       []),
                          rxn.Reactant(spec.DnaSpec('B', 1, spec.Domain(None, None, None)),
                                       [])],

                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       []),
                          rxn.Reactant(spec.DnaSpec('B', 1, spec.Domain(None, None, None)),
                                       []),
                          rxn.Reactant(spec.RnaSpec('B', 1, spec.Domain(None, None, None)),
                                       [])],

                         ),

        ReactionTestCase('A@0_trsl_BmRNA@1',
                         reaction_definitions['translation'],
                         rxn.Reaction(reaction_definitions['translation'],
                                      {'$x': rxn.specification_from_string('A@0'),
                                       '$y': rxn.specification_from_string('BmRNA@1')
                                       }),
                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       []),
                          rxn.Reactant(spec.RnaSpec('B', 1, spec.Domain(None, None, None)),
                                       [])],

                         [rxn.Reactant(spec.ProteinSpec('A', 0, spec.Domain(None, None, None)),
                                       []),
                          rxn.Reactant(spec.RnaSpec('B', 1, spec.Domain(None, None, None)),
                                       []),
                          rxn.Reactant(spec.ProteinSpec('B', 1, spec.Domain(None, None, None)),
                                       [])]
                         ),
            ]


def test_reactions(reaction_test_case):
    for the_case in reaction_test_case:
        is_reaction_correct(the_case)


def is_reaction_correct(the_case):
    actual_reaction = rxn.reaction_from_string(the_case.reaction_str)

    assert actual_reaction == the_case.expected_reaction and (actual_reaction == the_case.expected_reaction) is not None
    assert actual_reaction.definition == the_case.expected_definition and (actual_reaction.definition == the_case.expected_definition) is not None
    assert actual_reaction.reactants_lhs == the_case.expected_reactants_pre and (actual_reaction.reactants_lhs == the_case.expected_reactants_pre) is not None
    assert actual_reaction.reactants_rhs == the_case.expected_reactants_post and (actual_reaction.reactants_rhs == the_case.expected_reactants_post) is not None



def test_simple():

    reaction = rxn.reaction_from_string('A@0_[m]_bind_Agene@1_[n]')
    print(reaction)
    print(reaction.reactants_lhs)
    print(reaction.reactants_rhs)


#
