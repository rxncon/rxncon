#from rxncon.core.reaction import *
from collections import namedtuple
import pytest

import rxncon.core.reaction as rxn
import rxncon.core.specification as spec


DefinitionTestCase = namedtuple('DefinitionTestCase', ['definition', 'expected_name', 'expected_reactant_defs_pre',
                                                       'expected_reactant_defs_post', 'expected_reactants_defs',
                                                       'expected_representation_def', 'expected_variable_x',
                                                       'expected_variable_y'])


@pytest.fixture
def definition_test_case():
    return [
        DefinitionTestCase(rxn.ReactionDefinition('phosphorylation',
                                                  '$x_p+_$y',
                                                  {
                                                    '$x': (spec.ProteinSpecification, spec.SpecificationResolution.component),
                                                    '$y': (spec.Specification, spec.SpecificationResolution.residue)
                                                  },
                                                  '$x# + $y#$y-{0} -> $x# + $y#$y-{p}'
                                                  ),
                            'phosphorylation',
                            ['$x#', '$y#$y-{0}'],
                            ['$x#', '$y#$y-{p}'],
                            '$x# + $y#$y-{0} -> $x# + $y#$y-{p}',
                            '$x_p+_$y',
                            (spec.ProteinSpecification, spec.SpecificationResolution.component),
                            (spec.Specification, spec.SpecificationResolution.residue)),

        DefinitionTestCase(rxn.ReactionDefinition('phosphotransfer',
                                                  '$x_pt_$y',
                                                  {
                                                    '$x': (spec.Specification, spec.SpecificationResolution.residue),
                                                    '$y': (spec.Specification, spec.SpecificationResolution.residue)
                                                  },
                                                  '$x#$x-{p} + $y#$y-{0} -> $x#$x-{0} + $y#$y-{p}'
                                                  ),
                            'phosphotransfer',
                            ['$x#$x-{p}', '$y#$y-{0}'],
                            ['$x#$x-{0}', '$y#$y-{p}'],
                            '$x#$x-{p} + $y#$y-{0} -> $x#$x-{0} + $y#$y-{p}',
                            '$x_pt_$y',
                            (spec.Specification, spec.SpecificationResolution.residue),
                            (spec.Specification, spec.SpecificationResolution.residue)),

    DefinitionTestCase(rxn.ReactionDefinition('protein-protein-interaction',
                                              '$x_ppi_$y',
                                              {
                                               '$x': (spec.ProteinSpecification, spec.SpecificationResolution.domain),
                                                '$y': (spec.ProteinSpecification, spec.SpecificationResolution.domain)
                                              },
                                              '$x#$x--0 + $y#$y--0 <-> $x#$x--$y + $y#$x--$y'
                                             ),
                        'protein-protein-interaction',
                        ['$x#$x--0', '$y#$y--0'],
                        ['$x#$x--$y', '$y#$x--$y'],
                        '$x#$x--0 + $y#$y--0 <-> $x#$x--$y + $y#$x--$y',
                        '$x_ppi_$y',
                        (spec.ProteinSpecification, spec.SpecificationResolution.domain),
                        (spec.ProteinSpecification, spec.SpecificationResolution.domain)),

    DefinitionTestCase(rxn.ReactionDefinition('transcription',
                                              '$x_trsc_$y',
                                              {
                                                '$x': (spec.ProteinSpecification, spec.SpecificationResolution.component),
                                                '$y': (spec.ProteinSpecification, spec.SpecificationResolution.component)
                                              },
                                              '$x# + $y.gene# -> $x# + $y.gene# + $y.mRNA#'
                                             ),
                        'transcription',
                        ['$x#', '$y.gene#'],
                        ['$x#', '$y.gene#', '$y.mRNA#'],
                        '$x# + $y.gene# -> $x# + $y.gene# + $y.mRNA#',
                        '$x_trsc_$y',
                        (spec.ProteinSpecification, spec.SpecificationResolution.component),
                        (spec.ProteinSpecification, spec.SpecificationResolution.component)
    ),
    DefinitionTestCase(rxn.ReactionDefinition('translation',
                                              '$x_trsl_$y',
                                              {
                                                '$x': (spec.ProteinSpecification, spec.SpecificationResolution.component),
                                                '$y': (spec.ProteinSpecification, spec.SpecificationResolution.component)
                                              },
                                              '$x# + $y.mRNA# -> $x# + $y.mRNA# + $y#'
                                             ),
                        'translation',
                        ['$x#', '$y.mRNA#'],
                        ['$x#', '$y.mRNA#', '$y#'],
                        '$x# + $y.mRNA# -> $x# + $y.mRNA# + $y#',
                        '$x_trsl_$y',
                        (spec.ProteinSpecification, spec.SpecificationResolution.component),
                        (spec.ProteinSpecification, spec.SpecificationResolution.component)
                    ),
        DefinitionTestCase(rxn.ReactionDefinition('synthesis',
                                                  '$x_syn_$y',
                                                  {
                                                    '$x': (spec.ProteinSpecification, spec.SpecificationResolution.component),
                                                    '$y': (spec.ProteinSpecification, spec.SpecificationResolution.component)
                                                  },
                                                  '$x# -> $x# + $y#'
                                            ),
                       'synthesis',
                       ['$x#'],
                       ['$x#', '$y#'],
                       '$x# -> $x# + $y#',
                       '$x_syn_$y',
                       (spec.ProteinSpecification, spec.SpecificationResolution.component),
                       (spec.ProteinSpecification, spec.SpecificationResolution.component)
        ),

        DefinitionTestCase(rxn.ReactionDefinition('degradation',
                                                  '$x_deg_$y',
                                                  {
                                                      '$x': (spec.ProteinSpecification,spec.SpecificationResolution.component),
                                                      '$y': (spec.ProteinSpecification, spec.SpecificationResolution.component)
                                                  },
                                                  '$x# + $y# -> $x#'
                                                  ),
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



def test_simple():
    reaction = rxn.reaction_from_string(rxn.REACTION_DEFINITIONS, 'A_p+_B_[(y)]')
    print(reaction)
    print(reaction.reactants_pre)
    print(reaction.reactants_post)


