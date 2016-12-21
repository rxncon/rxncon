from enum import Enum
from typing import Tuple, List

from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, MolDef, Complex, \
    InitialCondition, Mol, Parameter, Rule, Observable


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


def bngl_from_rule_based_model(rule_based_model: RuleBasedModel, settings=BNGLSettings()) -> str:
    def header_str() -> str:
        return 'begin model'

    def molecule_types_str() -> str:
        molecule_types = [str_from_mol_def(mol_def) for mol_def in sorted(rule_based_model.mol_defs, key=lambda x: x.name)]
        return 'begin molecule types\n{0}\nend molecule types\n'.format('\n'.join(molecule_types))

    def seed_species_str() -> str:
        seeded_species = [str_from_initial_condition(initial_condition) for initial_condition in rule_based_model.initial_conditions]
        return 'begin seed species\n{0}\nend seed species\n'.format('\n'.join(sorted(seeded_species)))

    def parameters_str() -> str:
        parameters = [str_from_parameter(parameter) for parameter
                      in rule_based_model.parameters +
                      [ic.value for ic in rule_based_model.initial_conditions] +
                      rule_based_model.rate_parameters]
        return 'begin parameters\n{0}\nend parameters\n'.format('\n'.join(sorted(parameters)))

    def observables_str() -> str:
        observables = [str_from_observable(observable) for observable in rule_based_model.observables]
        return 'begin observables\n{0}\nend observables\n'.format('\n'.join(sorted(observables)))

    def reaction_rules_str() -> str:
        rules = [str_from_rule(rule) for rule in rule_based_model.rules]
        return 'begin reaction rules\n{0}\nend reaction rules\n'.format('\n'.join(rules))

    def footer_str() -> str:
        return 'end model\n\nsimulate_nf({{suffix=>nf,t_end=>100,n_steps=>50}});\n'\
            .format(settings.maximal_iteration,
                    settings.maximal_aggregate,
                    settings.simulation_method.value,
                    settings.simulation_time_end,
                    settings.simulation_time_steps)

    bngl_strs = [header_str(),
                 parameters_str(),
                 molecule_types_str(),
                 seed_species_str(),
                 observables_str(),
                 reaction_rules_str(),
                 footer_str()]

    return '\n'.join(bngl_str for bngl_str in bngl_strs if bngl_str)


def str_from_mol_def(mol_def: MolDef) -> str:
    def site_str(site_def: Tuple[str, List[str]]) -> str:
        return '~'.join([site_def[0]] + site_def[1])

    return '{0}({1})'.format(mol_def.name, ','.join(site_str(x) for x in mol_def.site_defs.items()))


def str_from_mol(mol: Mol) -> str:
    def site_str(site: str) -> str:
        s = site
        if site in mol.site_to_mod.keys():
            s += '~{}'.format(mol.site_to_mod[site])
        if site in mol.site_to_bond.keys():
            s += '!{}'.format(mol.site_to_bond[site]) if mol.site_to_bond[site] is not None else ''

        return s

    return '{0}({1})'.format(mol.name, ','.join(site_str(x) for x in mol.sites))


def str_from_complex(complex: Complex) -> str:
    return '.'.join(str_from_mol(mol) for mol in complex.mols)


def str_from_initial_condition(initial_condition: InitialCondition) -> str:
    value_str = initial_condition.value.name if initial_condition.value.name else initial_condition.value.value

    return '{0}\t{1}'.format(str_from_complex(initial_condition.complex), value_str)


def str_from_parameter(parameter: Parameter) -> str:
    assert parameter.name and parameter.value

    return '{0}\t\t{1}'.format(parameter.name, parameter.value)


def str_from_observable(observable: Observable) -> str:
    def clean(name: str) -> str:
        bad_chars = ['-', '[', ']']
        for bad_char in bad_chars:
            name = name.replace(bad_char, '')

        return name

    return 'Molecules\t{0}\t{1}'.format(clean(observable.name), str_from_complex(observable.complex))


def str_from_rule(rule: Rule) -> str:
    return '{0} -> {1}   {2}   # {3}'.format(' + '.join(str_from_complex(x) for x in rule.lhs),
                                             ' + '.join(str_from_complex(x) for x in rule.rhs),
                                             rule.rate.name if rule.rate.name else rule.rate.value,
                                             str(rule.parent_reaction))
