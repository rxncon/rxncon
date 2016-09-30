from enum import Enum
from typing import Tuple, List

from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, MolDef, Complex, SiteName, SiteModifier, \
    InitialCondition, Mol

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


def bngl_str_from_rule_based_model(rule_based_model: RuleBasedModel, settings=BNGLSettings()) -> str:
    def header_str() -> str:
        return 'begin model'

    def molecule_types_str() -> str:
        molecule_types = [str_from_mol_def(mol_def) for mol_def in sorted(rule_based_model.mol_defs)]
        return 'begin molecule types\n{0}\nend molecule types\n'.format('\n'.join(molecule_types))

    def seed_species_str() -> str:
        seeded_species = [str_from_initial_condition(initial_condition) for initial_condition in rule_based_model.initial_conditions]
        return 'begin seed species\n{0}\nend seed species\n'.format('\n'.join(seeded_species))

    def parameters_str() -> str:
        parameters = [string_from_parameter(parameter) for parameter in rule_based_model.parameters]
        return 'begin parameters\n{0}\nend parameters\n'.format('\n'.join(parameters))

    def observables_str() -> str:
        # @todo Add this.
        return ''

    def reaction_rules_str() -> str:
        rules = [string_from_rule(rule) for rule in rule_based_model.rules]
        return 'begin reaction rules\n{0}\nend reaction rules\n'.format('\n'.join(rules))

    def footer_str() -> str:
        return 'end model\n\ngenerate_network(max_iter=>{0}, max_agg=>{1})\nsimulate({{method=>\"{2}\",t_end=>{3},n_steps=>{4}}})\n'\
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
    bngl_strs = [bngl_str for bngl_str in bngl_strs if bngl_str]

    return '\n'.join(bngl_strs)


def str_from_mol_def(mol_def: MolDef) -> str:
    def site_str(site_def: Tuple[SiteName, List[SiteModifier]]) -> str:
        return '~'.join([site_def[0]] + site_def[1])

    return '{0}({1})'.format(mol_def.name, ','.join(site_str(x) for x in mol_def.site_defs.items()))


def str_from_mol(mol: Mol) -> str:
    def site_str(site: SiteName) -> str:
        site_str = site
        if mol.site_modifiers[site]:
            site_str += '~{}'.format(mol.site_modifiers[site])
        if mol.site_bonds[site]:



def str_from_complex(complex: Complex) -> str:
    return '.'.join(str_from_mol(mol) for mol in complex.mols)


def str_from_initial_condition(initial_condition: InitialCondition) -> str:
    value_str = initial_condition.value.name if initial_condition.value.name else initial_condition.value.value

    return '{0}\t{1}'.format(str_from_complex(initial_condition.complex), value_str)
