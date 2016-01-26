from typing import List, Optional, Tuple

import rxncon.simulation.rule_based.rule_based_model as rbm


class BNGLSystem:
    def __init__(self, rule_based_model: rbm.RuleBasedModel, maximal_iteration: Optional[int] = 1,
                 maximal_aggregate: Optional[int]=4, simulation_method: Optional[str] = "ode",
                 simulation_time_end: Optional[int] = 10, simulation_time_steps: Optional[int] = 100):
        self.rule_based_model = rule_based_model
        self.maximal_iteration = maximal_iteration  # maximum number of iterations of rule application
        self.maximal_aggregate = maximal_aggregate  # maximum number of molecules in one species
        # simulation method in BNGL we have
        # ode and ssa(kinetic Monte Carlo simulation using the Gillespie algorithm)
        self.simulation_method = simulation_method
        self.simulation_time_end = simulation_time_end
        self.simulation_time_steps = simulation_time_steps

    def to_string(self):
        bngl_system_strs = [self._header_string(),
                            self._parameters_string(),
                            self._molecule_types_string(),
                            self._seed_species_string(),
                            self._observables_string(),
                            self._reaction_rules_string(),
                            self._footer_string()]
        bngl_system_strs = [bngl_system_str for bngl_system_str in bngl_system_strs if bngl_system_str != ""]
        return '\n'.join(bngl_system_strs)

    def _header_string(self):
        return "begin model"

    def _molecule_types_string(self):
        molecule_types = [string_from_molecule_definition(molecule_definition) for molecule_definition in self.rule_based_model.molecule_definitions]
        return "begin molecule types\n{0}\nend molecule types\n".format("\n".join(molecule_types))

    def _seed_species_string(self):
        seeded_species = [string_from_initial_condition(inital_condition) for inital_condition in self.rule_based_model.initial_conditions]
        return "begin seed species\n{0}\nend seed species\n".format("\n".join(seeded_species))

    def _parameters_string(self):
        parameters = [string_from_parameter(parameter) for parameter in self.rule_based_model.parameters]
        return "begin parameters\n{0}\nend parameters\n".format("\n".join(parameters))

    def _observables_string(self):
        return ""

    def _reaction_rules_string(self):
        rules = [string_from_rule(rule) for rule in self.rule_based_model.rules]
        return "begin reaction rules\n{0}\nend reaction rules\n".format("\n".join(rules))

    def _footer_string(self):
        """
        The generate_network command directs BioNetGen to generate a network of species and reactions through iterative
        application of the rules starting from the set of seed species. At each step in this iterative process, rules
        are applied to the existing set of chemical species to generate new reactions. Following rule application,
        the species appearing as products in the new reactions are checked to determine whether they correspond to
        existing species in the network. If no new species are found, network generation terminates.

        max_iter: upper limit on the number of iterations of rule application
        max_agg: upper limit on the number molecules in an aggregate
        """

        return 'end model\n\ngenerate_network(max_iter=>{0}, max_agg=>{1})\n' \
               'simulate({{method=>"{2}",t_end=>{3},n_steps=>{4}}})'.format(self.maximal_iteration, self.maximal_aggregate,
                                                                         self.simulation_method, self.simulation_time_end,
                                                                         self.simulation_time_steps)


# MOLECULE DEF / SPEC
def string_from_molecule_definition(molecule_definition: rbm.MoleculeDefinition) -> str:
    if not molecule_definition.modification_definitions and not molecule_definition.association_definitions and not\
            molecule_definition.localization_definitions:
        return molecule_definition.name

    domain_strs = [','.join([string_from_localization_definition(x) for x in molecule_definition.localization_definitions]),
                   ','.join([string_from_modification_definition(x) for x in molecule_definition.modification_definitions]),
                   ','.join([string_from_association_definition(x) for x in molecule_definition.association_definitions])]
    domain_strs = [domain_str for domain_str in domain_strs if domain_str != ""]
    return "{0}({1})".format(molecule_definition.name,",".join(domain_strs))


def string_from_molecule_specification(molecule_specification: rbm.MoleculeSpecification, rxn_binding: Optional[Tuple[int, rbm.AssociationSpecification]] = None) -> str:
    if not molecule_specification.modification_specifications and not molecule_specification.association_specifications and not\
            molecule_specification.localization_specifications:
        return molecule_specification.molecule_definition.name

    domain_strs = [','.join(string_from_localization_specification(x) for x in molecule_specification.localization_specifications),
                   ','.join(string_from_modification_specification(x) for x in molecule_specification.modification_specifications),
                   ','.join(string_from_association_specification(x) if rxn_binding is None else string_from_bound_associaton_specification(x, rxn_binding) for x in molecule_specification.association_specifications)]

    domain_strs = [domain_str for domain_str in domain_strs if domain_str != ""]
    return "{0}({1})".format(molecule_specification.molecule_definition.name, ",".join(domain_strs))


def string_from_bound_associaton_specification(assocDom: rbm.AssociationSpecification, rxn_binding: Tuple[int, rbm.AssociationSpecification]):
    if assocDom == rxn_binding[1]:
        return "{0}!{1}".format(string_from_association_specification(assocDom), rxn_binding[0])
    return assocDom


# MODIFICATION DEF / SPEC
def string_from_modification_definition(modification_definition: rbm.ModificationDefinition) -> str:
    return '{0}~{1}'.format(modification_definition.domain_name, '~'.join(modification_definition.valid_modifiers))


def string_from_modification_specification(modification_specification: rbm.ModificationSpecification) -> str:
    return "{0}~{1}".format(modification_specification.modification_definition.domain_name, modification_specification.value)


# ASSOCIATION DEF / SPEC
def string_from_association_definition(association_definition: rbm.AssociationDefinition) -> str:
    return association_definition.domain_name


def string_from_association_specification(association_specification: rbm.AssociationSpecification) -> str:
    return association_specification.association_definition.domain_name


# LOCALIZATION DEF / SPEC
def string_from_localization_definition(localization_definition: rbm.LocalizationDefinition) -> str:
    return 'loc~{0}'.format('~'.join(localization_definition.compartments))


def string_from_localization_specification(localization_specification: rbm.LocalizationSpecification) -> str:
    return 'loc~' + localization_specification.current_compartment


# REACTANTS
def string_from_molecule_reactant(molecule_reactant: rbm.MoleculeReactant) -> str:
    return string_from_molecule_specification(molecule_reactant.molecule_specification)


def string_from_complex_reactant(complex_reactant: rbm.ComplexReactant) -> str:
    def __is_association_domain_bound(i: int, complex_part: rbm.MoleculeSpecification, rxn_number: int, binding: rbm.Binding):
        if binding.left_partner[0] == i:
            return string_from_molecule_specification(complex_part,(rxn_number,binding.left_partner[1]))
        elif binding.right_partner[0] == i:
             return string_from_molecule_specification(complex_part,(rxn_number,binding.right_partner[1]))

    complex = [__is_association_domain_bound(i, complex_part, rxn_number, binding)
               for i, complex_part in enumerate(complex_reactant.complex_parts)
               for rxn_number, binding in enumerate(complex_reactant.complex_bindings)]
    return ".".join(complex)



def string_from_rule(rule: rbm.Rule):
    def __helper(reactant: rbm.Reactant):
        if isinstance(reactant, rbm.MoleculeReactant):
            return string_from_molecule_reactant(reactant)
        elif isinstance(reactant, rbm.ComplexReactant):
            return string_from_complex_reactant(reactant)
        else:
            raise "Reactant is not defined"

    left_hand_side = [__helper(reactant) for reactant in rule.left_hand_side]
    right_hand_side = [__helper(reactant) for reactant in rule.right_hand_side]
    kinetic_parameters = [parameter.name for parameter in rule.parameters]
    return "{0} {1} {2} {3}".format(" + ".join(left_hand_side), rule.arrow_type.value, " + ".join(right_hand_side), ",".join(kinetic_parameters))


def string_from_parameter(parameter: rbm.Parameter):
    return "{0} {1}".format(parameter.name, parameter.value)


def string_from_initial_condition(initial_condition: rbm.InitialCondition):
    return  "{0} {1}".format(string_from_molecule_specification(initial_condition.molecule_specification), initial_condition.value)