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
    return {'phosphorylation': rxn.ReactionDefinition('phosphorylation',
                                                      '$x_p+_$y',
                                                      {
                                                          '$x': (spec.ProteinSpecification, spec.SpecificationResolution.component),
                                                          '$y': (spec.Specification, spec.SpecificationResolution.residue)
                                                      },
                                                      '$x# + $y#$y-{0} -> $x# + $y#$y-{p}'
                                                      ),

            'phosphotransfer': rxn.ReactionDefinition('phosphotransfer',
                                                      '$x_pt_$y',
                                                      {
                                                        '$x': (spec.Specification, spec.SpecificationResolution.residue),
                                                        '$y': (spec.Specification, spec.SpecificationResolution.residue)
                                                      },
                                                      '$x#$x-{p} + $y#$y-{0} -> $x#$x-{0} + $y#$y-{p}'
                                                      ),
            'protein-protein-interaction': rxn.ReactionDefinition('protein-protein-interaction',
                                                                  '$x_ppi_$y',
                                                                  {
                                                                    '$x': (spec.ProteinSpecification, spec.SpecificationResolution.domain),
                                                                   '$y': (spec.ProteinSpecification, spec.SpecificationResolution.domain)
                                                                  },
                                                                  '$x#$x--0 + $y#$y--0 <-> $x#$x--$y + $y#$x--$y'
                                                                 ),

            'intra-protein-interaction': rxn.ReactionDefinition('intra-protein-interaction',
                                                            '$x_ipi_$y',
                                                            {
                                                                '$x': (spec.ProteinSpecification, spec.SpecificationResolution.domain),
                                                                '$y': (spec.ProteinSpecification, spec.SpecificationResolution.domain)
                                                            },
                                                            '$x#$x--0,$y--0 -> $x#$x--$y.domain'
                                                            ),

            'transcription': rxn.ReactionDefinition('transcription',
                                                   '$x_trsc_$y',
                                                    {
                                                        '$x': (spec.ProteinSpecification, spec.SpecificationResolution.component),
                                                        '$y': (spec.Specification, spec.SpecificationResolution.component)
                                                    },
                                                   '$x# + $y.gene# -> $x# + $y.gene# + $y.mRNA#'
                                                   ),

            'translation': rxn.ReactionDefinition('translation',
                                                  '$x_trsl_$y',
                                                  {
                                                    '$x': (spec.ProteinSpecification, spec.SpecificationResolution.component),
                                                    '$y': (spec.Specification, spec.SpecificationResolution.component)
                                                  },
                                                  '$x# + $y.mRNA# -> $x# + $y.mRNA# + $y#'
                                                 ),

            'synthesis': rxn.ReactionDefinition('synthesis',
                                                '$x_syn_$y',
                                                {
                                                  '$x': (spec.ProteinSpecification, spec.SpecificationResolution.component),
                                                  '$y': (spec.ProteinSpecification, spec.SpecificationResolution.component)
                                                },
                                                '$x# -> $x# + $y#'
                                                ),

            'degradation': rxn.ReactionDefinition('degradation',
                                                  '$x_deg_$y',
                                                  {
                                                      '$x': (spec.ProteinSpecification,spec.SpecificationResolution.component),
                                                      '$y': (spec.ProteinSpecification, spec.SpecificationResolution.component)
                                                  },
                                                  '$x# + $y# -> $x#'
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
                           (spec.ProteinSpecification, spec.SpecificationResolution.component),
                           (spec.Specification, spec.SpecificationResolution.residue)),

        DefinitionTestCase(reaction_definitions['phosphotransfer'],
                            'phosphotransfer',
                           ['$x#$x-{p}', '$y#$y-{0}'],
                           ['$x#$x-{0}', '$y#$y-{p}'],
                            '$x#$x-{p} + $y#$y-{0} -> $x#$x-{0} + $y#$y-{p}',
                            '$x_pt_$y',
                           (spec.Specification, spec.SpecificationResolution.residue),
                           (spec.Specification, spec.SpecificationResolution.residue)),

    DefinitionTestCase(reaction_definitions['protein-protein-interaction'],
                        'protein-protein-interaction',
                       ['$x#$x--0', '$y#$y--0'],
                       ['$x#$x--$y', '$y#$x--$y'],
                        '$x#$x--0 + $y#$y--0 <-> $x#$x--$y + $y#$x--$y',
                        '$x_ppi_$y',
                       (spec.ProteinSpecification, spec.SpecificationResolution.domain),
                       (spec.ProteinSpecification, spec.SpecificationResolution.domain)),

    DefinitionTestCase(reaction_definitions['transcription'],
                        'transcription',
                       ['$x#', '$y.gene#'],
                       ['$x#', '$y.gene#', '$y.mRNA#'],
                        '$x# + $y.gene# -> $x# + $y.gene# + $y.mRNA#',
                        '$x_trsc_$y',
                       (spec.ProteinSpecification, spec.SpecificationResolution.component),
                       (spec.Specification, spec.SpecificationResolution.component)
                       ),
    DefinitionTestCase(reaction_definitions['translation'],
                        'translation',
                       ['$x#', '$y.mRNA#'],
                       ['$x#', '$y.mRNA#', '$y#'],
                        '$x# + $y.mRNA# -> $x# + $y.mRNA# + $y#',
                        '$x_trsl_$y',
                       (spec.ProteinSpecification, spec.SpecificationResolution.component),
                       (spec.Specification, spec.SpecificationResolution.component)
                       ),
        DefinitionTestCase(reaction_definitions['synthesis'],
                       'synthesis',
                           ['$x#'],
                           ['$x#', '$y#'],
                       '$x# -> $x# + $y#',
                       '$x_syn_$y',
                           (spec.ProteinSpecification, spec.SpecificationResolution.component),
                           (spec.ProteinSpecification, spec.SpecificationResolution.component)
                           ),

        DefinitionTestCase(reaction_definitions['degradation'],
                           'degradation',
                           ['$x#', '$y#'],
                           ['$x#'],
                           '$x# + $y# -> $x#',
                           '$x_deg_$y',
                           (spec.ProteinSpecification, spec.SpecificationResolution.component),
                           (spec.ProteinSpecification, spec.SpecificationResolution.component)
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
        ReactionTestCase('A_p+_B_[(y)]',
                         reaction_definitions['phosphorylation'],
                         rxn.Reaction(reaction_definitions['phosphorylation'],
                                      {'$x': rxn.specification_from_string('A'),
                                       '$y': rxn.specification_from_string('B_[(y)]')
                                       }),

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('B_[(y)]-{0}')])],

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('B_[(y)]-{P}')])],

                         ),

        ReactionTestCase('A_[(x)]_pt_B_[(y)]',  # todo: residue has to be defined is this correct?
                         reaction_definitions['phosphotransfer'],
                         rxn.Reaction(reaction_definitions['phosphotransfer'],
                                      {'$x': rxn.specification_from_string('A_[(x)]'),
                                       '$y': rxn.specification_from_string('B_[(y)]')
                                       }),
                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('A_[(x)]-{P}')]),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('B_[(y)]-{0}')])],

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('A_[(x)]-{0}')]),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('B_[(y)]-{P}')])],
                         ),

        ReactionTestCase('A_[b]_ppi_B_[a]',
                         reaction_definitions['protein-protein-interaction'],
                         rxn.Reaction(reaction_definitions['protein-protein-interaction'],
                                      {'$x': rxn.specification_from_string('A_[b]'),
                                       '$y': rxn.specification_from_string('B_[a]')
                                       }),
                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('A_[b]--0')]),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('B_[a]--0')])],

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('A_[b]--B_[a]')]),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('A_[b]--B_[a]')])],

                         ),
        ReactionTestCase('A_[n]_ipi_A_[m]',
                         reaction_definitions['intra-protein-interaction'],
                         rxn.Reaction(reaction_definitions['intra-protein-interaction'],
                                      {'$x': rxn.specification_from_string('A_[n]'),
                                       '$y': rxn.specification_from_string('A_[m]')
                                       }),
                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('A_[n]--0'), sta.state_from_string('A_[m]--0')])],

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       [sta.state_from_string('A_[n]--[m]')])],

                         ),

            # ReactionTestCase('A_[b]_ppi_B_[d/s]',
            #                  reaction_definitions['protein-protein-interaction'],
            #                  rxn.Reaction(reaction_definitions['protein-protein-interaction'],
            #                               { '$x': rxn.specification_from_string('A_[b]'),
            #                                 '$y': rxn.specification_from_string('B_[d/s]')
            #                                 }),
            #                  [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
            #                                [sta.state_from_string('A_[b]--0')]),
            #                   rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
            #                                [sta.state_from_string('B_[d/s]--0')])],
            #
            #                  [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
            #                                [sta.state_from_string('A_[b]--B_[a]')]),
            #                   rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
            #                                [sta.state_from_string('A_[d/s]--B_[a]')])],
            #
            #                  ),

        ReactionTestCase('A_trsc_B',
                         reaction_definitions['transcription'],
                         rxn.Reaction(reaction_definitions['transcription'],
                                      {'$x': rxn.specification_from_string('A'),
                                       '$y': rxn.specification_from_string('B')
                                       }),
                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.DnaSpecification('B', spec.DomainResolution(None, None, None)),
                                       [])],

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.DnaSpecification('B', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.RnaSpecification('B', spec.DomainResolution(None, None, None)),
                                       [])],

                         ),

        ReactionTestCase('A_trsl_B',
                         reaction_definitions['translation'],
                         rxn.Reaction(reaction_definitions['translation'],
                                      {'$x': rxn.specification_from_string('A'),
                                       '$y': rxn.specification_from_string('B')
                                       }),
                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.RnaSpecification('B', spec.DomainResolution(None, None, None)),
                                       [])],

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.RnaSpecification('B', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       [])]
                         ),
            #
        ReactionTestCase('A_syn_B',
                         reaction_definitions['synthesis'],
                         rxn.Reaction(reaction_definitions['synthesis'],
                                      {'$x': rxn.specification_from_string('A'),
                                       '$y': rxn.specification_from_string('B')
                                       }),
                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       [])],

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       []),
                          ]
                         ),
        #
        ReactionTestCase('A_deg_B',
                         reaction_definitions['degradation'],
                         rxn.Reaction(reaction_definitions['degradation'],
                                      {'$x': rxn.specification_from_string('A'),
                                       '$y': rxn.specification_from_string('B')
                                       }),
                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       []),
                          rxn.Reactant(spec.ProteinSpecification('B', spec.DomainResolution(None, None, None)),
                                       [])],

                         [rxn.Reactant(spec.ProteinSpecification('A', spec.DomainResolution(None, None, None)),
                                       [])
                          ]
                         )
            ]


def test_reactions(reaction_test_case):
    for the_case in reaction_test_case:
        is_reaction_correct(the_case)


def is_reaction_correct(the_case):
    actual_reaction = rxn.reaction_from_string(the_case.reaction_str)

    assert actual_reaction == the_case.expected_reaction
    assert actual_reaction.definition == the_case.expected_definition
    assert actual_reaction.reactants_pre == the_case.expected_reactants_pre
    assert actual_reaction.reactants_post == the_case.expected_reactants_post


#
def test_simple():

    reaction = rxn.reaction_from_string('A_[m]_bind_A_[n]')
    print(reaction)
    print(reaction.reactants_pre)
    print(reaction.reactants_post)



