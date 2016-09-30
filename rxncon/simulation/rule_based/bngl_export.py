from enum import Enum
from typing import Optional

from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel

class BNGLSimulationMethods(Enum):
    ODE = 'ode'
    SSA = 'ssa'


class BNGLSettings:
    def __init__(self):
        self.maximal_iteration = 1
        self.maximal_aggregate = 4
        self.simulation_method = BNGLSimulationMethods.ODE
        self.simulation_time_end = 10
        self.simulation_time_steps = 100


def bngl_str_from_rule_based_model(rule_based_model: RuleBasedModel, setting: Optional[BNGLSettings]):
    def __init__(self, rule_based_model: RuleBasedModel, settings: BNGLSettings):
        self.rule_based_model = rule_based_model
        self.settings = settings

    def to_string(self) -> str:
        bngl_system_strs = [self._header_string(),
                            self._parameters_string(),
                            self._molecule_types_string(),
                            self._seed_species_string(),
                            self._observables_string(),
                            self._reaction_rules_string(),
                            self._footer_string()]
        bngl_system_strs = [bngl_system_str for bngl_system_str in bngl_system_strs if bngl_system_str]

        return '\n'.join(bngl_system_strs)

    def _header_string(self) -> str:
        return 'begin model'

    def _molecule_types_string(self) -> str:
        molecule_types = [string_from_mol_def(molecule_definition) for molecule_definition in sorted(self.rule_based_model.mol_defs)]
        return 'begin molecule types\n{0}\nend molecule types\n'.format('\n'.join(molecule_types))

    def _seed_species_string(self) -> str:
        seeded_species = [string_from_initial_condition(initial_condition) for initial_condition in self.rule_based_model.initial_conditions]
        return 'begin seed species\n{0}\nend seed species\n'.format('\n'.join(seeded_species))

    def _parameters_string(self) -> str:
        parameters = [string_from_parameter(parameter) for parameter in self.rule_based_model.parameters]
        return 'begin parameters\n{0}\nend parameters\n'.format('\n'.join(parameters))

    def _observables_string(self) -> str:
        # @todo Add this.
        return ''

    def _reaction_rules_string(self) -> str:
        rules = [string_from_rule(rule) for rule in self.rule_based_model.rules]
        return 'begin reaction rules\n{0}\nend reaction rules\n'.format('\n'.join(rules))

    def _footer_string(self) -> str:
        return 'end model\n\ngenerate_network(max_iter=>{0}, max_agg=>{1})\nsimulate({{method=>\"{2}\",t_end=>{3},n_steps=>{4}}})'\
            .format(self.settings.maximal_iteration,
                    self.settings.maximal_aggregate,
                    self.settings.simulation_method.value,
                    self.settings.simulation_time_end,
                    self.settings.simulation_time_steps)


def string_from_mol_def(molecule_definition: MoleculeDefinition) -> str:
    if not molecule_definition.mod_defs and not molecule_definition.ass_defs and not\
            molecule_definition.loc_def:
        return molecule_definition.spec

    definition_strings = [string_from_loc_prop_def(sorted(molecule_definition.loc_def))
                          if molecule_definition.loc_def else None,
                          ','.join([string_from_mod_prop_def(x) for x in sorted(molecule_definition.mod_defs)]),
                          ','.join([string_from_ass_prop_def(x) for x in sorted(molecule_definition.ass_defs)])]

    return '{0}({1})'.format(molecule_definition.spec, ','.join(x for x in definition_strings if x))


def string_from_molecule(molecule_specification: Molecule) -> str:
    if not molecule_specification.mod_props and not molecule_specification.ass_props and not\
            molecule_specification.loc_prop:
        return molecule_specification.mol_def.name

    specification_strings = [string_from_loc_prop(sorted(molecule_specification.loc_prop))
                             if molecule_specification.loc_prop else None,
                             ','.join(string_from_mod_prop(x) for x in sorted(molecule_specification.mod_props)),
                             ','.join(string_from_ass_prop(x) for x in sorted(molecule_specification.ass_props))]

    return '{0}({1})'.format(molecule_specification.mol_def.name, ','.join(x for x in specification_strings if x))


def string_from_mod_prop_def(modification_definition: ModificationPropertyDefinition) -> str:
    return '{0}~{1}'.format(modification_definition.spec, '~'.join(sorted(modification_definition.valid_modifiers)))


def string_from_mod_prop(modification_specification: ModificationProperty) -> str:
    return '{0}~{1}'.format(modification_specification.prop_def.domain, modification_specification.modifier)


def string_from_ass_prop_def(association_definition: AssociationPropertyDefinition) -> str:
    return association_definition.spec


def string_from_ass_prop(association_specification: AssociationProperty) -> str:
    if association_specification.occupation_status == OccupationStatus.not_specified:
        return association_specification.property_def.domain + '!?'

    elif association_specification.occupation_status == OccupationStatus.occupied_unknown_partner:
        return association_specification.property_def.domain + '!+'

    elif association_specification.occupation_status == OccupationStatus.not_occupied:
        return association_specification.property_def.domain

    elif association_specification.occupation_status == OccupationStatus.occupied_known_partner:
        raise NotImplementedError('AssociationSpecification with domain occupied by known partner cannot be stringified outside complex.')

    else:
        raise AssertionError('Unknown occupation status for association spec {0}'.format(association_specification))


def string_from_loc_prop_def(localization_definition: LocalizationPropertyDefinition) -> str:
    return 'loc~{0}'.format('~'.join(localization_definition.valid_compartments))


def string_from_loc_prop(localization_specification: LocalizationProperty) -> str:
    return 'loc~' + localization_specification.compartment


def string_from_complex(reactant: Reactant) -> str:
    if isinstance(reactant, MoleculeReactant):
        return string_from_molecule_reactant(reactant)

    elif isinstance(reactant, ComplexReactant):
        return string_from_complex_reactant(sorted(reactant))

    else:
        raise ValueError('Reactant {0} is neither MoleculeReactant nor ComplexReactant'.format(str(reactant)))


def string_from_rule(rule: rxncon.semantics.rule.Rule) -> str:
    left_hand_side = [string_from_complex(reactant) for reactant in sorted(rule.lhs)]
    right_hand_side = [string_from_complex(reactant) for reactant in sorted(rule.rhs)]
    kinetic_parameters = [parameter.name for parameter in sorted(rule.rates)]
    return '{0} {1} {2} {3}'.format(' + '.join(left_hand_side), rule.arrow_type.value, ' + '.join(right_hand_side), ','.join(kinetic_parameters))


def string_from_parameter(parameter: Parameter) -> str:
    return '{0} {1}'.format(parameter.name, parameter.value)


def string_from_initial_condition(initial_condition: InitialCondition) -> str:
    return '{0} {1}'.format(string_from_molecule(initial_condition.molecule), initial_condition.value)
