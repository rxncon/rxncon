from typing import List, Optional, Tuple
from collections import defaultdict
from enum import Enum

import rxncon.simulation.rule_based.rule_based_model as rbm


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


class BNGLSystem:
    def __init__(self, rule_based_model: rbm.RuleBasedModel, settings: BNGLSettings):
        self.rule_based_model = rule_based_model
        self.settings = settings

    def to_string(self):
        bngl_system_strs = [self._header_string(),
                            self._parameters_string(),
                            self._molecule_types_string(),
                            self._seed_species_string(),
                            self._observables_string(),
                            self._reaction_rules_string(),
                            self._footer_string()]
        bngl_system_strs = [bngl_system_str for bngl_system_str in bngl_system_strs if bngl_system_str != '']

        return '\n'.join(bngl_system_strs)

    def _header_string(self):
        return 'begin model'

    def _molecule_types_string(self):
        molecule_types = [string_from_molecule_definition(molecule_definition) for molecule_definition in self.rule_based_model.molecule_defs]
        return 'begin molecule types\n{0}\nend molecule types\n'.format('\n'.join(molecule_types))

    def _seed_species_string(self):
        seeded_species = [string_from_initial_condition(initial_condition) for initial_condition in self.rule_based_model.initial_conditions]
        return 'begin seed species\n{0}\nend seed species\n'.format('\n'.join(seeded_species))

    def _parameters_string(self):
        parameters = [string_from_parameter(parameter) for parameter in self.rule_based_model.parameters]
        return 'begin parameters\n{0}\nend parameters\n'.format('\n'.join(parameters))

    def _observables_string(self):
        # @todo Add this.
        return ''

    def _reaction_rules_string(self):
        rules = [string_from_rule(rule) for rule in self.rule_based_model.rules]
        return 'begin reaction rules\n{0}\nend reaction rules\n'.format('\n'.join(rules))

    def _footer_string(self):
        return 'end model\n\ngenerate_network(max_iter=>{0}, max_agg=>{1})\nsimulate({{method=>\"{2}\",t_end=>{3},n_steps=>{4}}})'\
            .format(self.settings.maximal_iteration,
                    self.settings.maximal_aggregate,
                    self.settings.simulation_method.value,
                    self.settings.simulation_time_end,
                    self.settings.simulation_time_steps)


# MOLECULE DEF / SPEC
def string_from_molecule_definition(molecule_definition: rbm.MoleculeDefinition) -> str:
    if not molecule_definition.modification_defs and not molecule_definition.association_defs and not\
            molecule_definition.localization_def:
        return molecule_definition.name

    definition_strings = [string_from_localization_definition(molecule_definition.localization_def)
                          if molecule_definition.localization_def else None,
                          ','.join([string_from_modification_definition(x) for x in molecule_definition.modification_defs]),
                          ','.join([string_from_association_definition(x) for x in molecule_definition.association_defs])]

    return '{0}({1})'.format(molecule_definition.name, ','.join(x for x in definition_strings if x))


def string_from_molecule_specification(molecule_specification: rbm.MoleculeSpecification) -> str:
    if not molecule_specification.modification_specs and not molecule_specification.association_specs and not\
            molecule_specification.localization_spec:
        return molecule_specification.molecule_def.name

    specification_strings = [string_from_localization_specification(molecule_specification.localization_spec)
                             if molecule_specification.localization_spec else None,
                             ','.join(string_from_modification_specification(x) for x in molecule_specification.modification_specs),
                             ','.join(string_from_association_specification(x) for x in molecule_specification.association_specs)]

    return '{0}({1})'.format(molecule_specification.molecule_def.name, ','.join(x for x in specification_strings if x))


# MODIFICATION DEF / SPEC
def string_from_modification_definition(modification_definition: rbm.ModificationDefinition) -> str:
    return '{0}~{1}'.format(modification_definition.domain, '~'.join(modification_definition.valid_modifiers))


def string_from_modification_specification(modification_specification: rbm.ModificationSpecification) -> str:
    return '{0}~{1}'.format(modification_specification.modification_def.domain, modification_specification.modifier)


# ASSOCIATION DEF / SPEC
def string_from_association_definition(association_definition: rbm.AssociationDefinition) -> str:
    return association_definition.domain


def string_from_association_specification(association_specification: rbm.AssociationSpecification) -> str:
    if association_specification.occupation_status == rbm.OccupationStatus.not_specified:
        return association_specification.association_def.domain + '!?'

    elif association_specification.occupation_status == rbm.OccupationStatus.occupied_unknown_partner:
        return association_specification.association_def.domain + '!+'

    elif association_specification.occupation_status == rbm.OccupationStatus.unoccupied:
        return association_specification.association_def.domain

    elif association_specification.occupation_status == rbm.OccupationStatus.occupied_known_partner:
        raise NotImplementedError('AssociationSpecification with domain occupied by known partner cannot be stringified outside complex.')

    else:
        raise AssertionError

# LOCALIZATION DEF / SPEC
def string_from_localization_definition(localization_definition: rbm.LocalizationDefinition) -> str:
    return 'loc~{0}'.format('~'.join(localization_definition.valid_compartments))


def string_from_localization_specification(localization_specification: rbm.LocalizationSpecification) -> str:
    return 'loc~' + localization_specification.compartment


# REACTANTS
def string_from_reactant(reactant: rbm.Reactant) -> str:
    if isinstance(reactant, rbm.MoleculeReactant):
        return string_from_molecule_reactant(reactant)

    elif isinstance(reactant, rbm.ComplexReactant):
        return string_from_complex_reactant(reactant)

    else:
        raise ValueError('Reactant {0} is neither MoleculeReactant nor ComplexReactant'.format(str(reactant)))


def string_from_molecule_reactant(molecule_reactant: rbm.MoleculeReactant) -> str:
    return string_from_molecule_specification(molecule_reactant.molecule_specification)


def string_from_complex_reactant(complex_reactant: rbm.ComplexReactant) -> str:
    index_to_assoc_to_binding_num = defaultdict(dict)

    for i, binding in enumerate(complex_reactant.bindings):
        for partner in (binding.left_partner, binding.right_partner):
            index = partner[0]
            assoc = partner[1]
            index_to_assoc_to_binding_num[index][assoc] = i

    complex_strings = []

    for i, molecule in enumerate(complex_reactant.molecules):
        bound_assocs = index_to_assoc_to_binding_num[i].keys()

        assoc_strings = []
        for assoc_spec in molecule.association_specs:
            if assoc_spec in bound_assocs:
                assoc_strings.append('{0}!{1}'.format(assoc_spec.association_def.domain,
                                                      index_to_assoc_to_binding_num[i][assoc_spec]))
            else:
                assoc_strings.append(string_from_association_specification(assoc_spec))

        specification_strings = [string_from_localization_specification(molecule.localization_spec)
                                 if molecule.localization_spec else None,
                                 ','.join(string_from_modification_specification(x) for x in molecule.modification_specs),
                                 ','.join(assoc_strings)]

        complex_strings.append('{0}({1})'.format(molecule.molecule_def.name, ','.join(x for x in specification_strings if x)))

    return '.'.join(complex_strings)


def string_from_rule(rule: rbm.Rule):
    left_hand_side = [string_from_reactant(reactant) for reactant in rule.left_hand_side]
    right_hand_side = [string_from_reactant(reactant) for reactant in rule.right_hand_side]
    kinetic_parameters = [parameter.name for parameter in rule.rates]
    return '{0} {1} {2} {3}'.format(' + '.join(left_hand_side), rule.arrow_type.value, ' + '.join(right_hand_side), ','.join(kinetic_parameters))


def string_from_parameter(parameter: rbm.Parameter):
    return '{0} {1}'.format(parameter.name, parameter.value)


def string_from_initial_condition(initial_condition: rbm.InitialCondition):
    return '{0} {1}'.format(string_from_molecule_specification(initial_condition.molecule_specification), initial_condition.value)