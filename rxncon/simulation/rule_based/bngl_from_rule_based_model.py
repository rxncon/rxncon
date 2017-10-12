"""Module containing the function bngl_from_rule_based_model, which translates a rule-based model
to BNGL."""

from typing import List

from rxncon.simulation.rule_based.rule_based_model import RuleBasedModel, MolDef, Complex, \
    InitialCondition, Mol, Parameter, Rule, Observable


def bngl_from_rule_based_model(rule_based_model: RuleBasedModel) -> str:
    def header_str() -> str:
        return 'begin model'

    def molecule_types_str() -> str:
        molecule_types = [_str_from_mol_def(mol_def) for mol_def in sorted(rule_based_model.mol_defs, key=lambda x: x.name)]
        return 'begin molecule types\n{0}\nend molecule types\n'.format('\n'.join(molecule_types))

    def seed_species_str() -> str:
        seeded_species = [_str_from_initial_condition(initial_condition) for initial_condition in rule_based_model.initial_conditions]
        return 'begin seed species\n{0}\nend seed species\n'.format('\n'.join(sorted(seeded_species)))

    def parameters_str() -> str:
        parameters = [_str_from_parameter(parameter) for parameter
                      in rule_based_model.parameters +
                      [ic.value for ic in rule_based_model.initial_conditions] +
                      rule_based_model.rate_parameters]
        return 'begin parameters\n{0}\nend parameters\n'.format('\n'.join(sorted(parameters)))

    def observables_str() -> str:
        observables = [_str_from_observable(observable) for observable in rule_based_model.observables]
        return 'begin observables\n{0}\nend observables\n'.format('\n'.join(sorted(observables)))

    def reaction_rules_str() -> str:
        rules = [_str_from_rule(rule, i) for i, rule in enumerate(rule_based_model.rules)]
        return 'begin reaction rules\n{0}\nend reaction rules\n'.format('\n'.join(rules))

    def footer_str() -> str:
        return 'end model\n\nsimulate_nf({t_end=>100,n_steps=>10});\n'

    bngl_strs = [header_str(),
                 parameters_str(),
                 molecule_types_str(),
                 seed_species_str(),
                 observables_str(),
                 reaction_rules_str(),
                 footer_str()]

    return '\n'.join(bngl_str for bngl_str in bngl_strs if bngl_str)


def _str_from_mol_def(mol_def: MolDef) -> str:
    def site_str(site_name: str, site_def: List[str]) -> str:
        return '~'.join([site_name] + sorted(site_def))

    return '{0}({1})'.format(mol_def.name, ','.join(site_str(site, mol_def.site_defs[site])
                                                    for site in sorted(mol_def.site_defs)))


def _str_from_mol(mol: Mol) -> str:
    def full_site_str(site: str) -> str:
        site_str = site
        if site in mol.site_to_mod.keys():
            site_str += '~{}'.format(mol.site_to_mod[site])
        if site in mol.site_to_bond.keys():
            site_str += '!{}'.format(mol.site_to_bond[site]) if mol.site_to_bond[site] is not None else ''

        return site_str

    return '{0}({1})'.format(mol.name, ','.join(full_site_str(x) for x in sorted(mol.sites)))


def _str_from_complex(complex: Complex) -> str:  # pylint: disable=redefined-builtin
    return '.'.join(_str_from_mol(mol) for mol in complex.mols)


def _str_from_initial_condition(initial_condition: InitialCondition) -> str:
    value_str = initial_condition.value.name if initial_condition.value.name else initial_condition.value.value

    return '{0}\t{1}'.format(_str_from_complex(initial_condition.complex), value_str)


def _str_from_parameter(parameter: Parameter) -> str:
    assert parameter.name and parameter.value
    if parameter.description:
        return '{0:<10}{1}\t\t#  {2}'.format(parameter.name, parameter.value, parameter.description)
    else:
        return '{0:<10}{1}'.format(parameter.name, parameter.value)


def _str_from_observable(observable: Observable) -> str:
    def clean(name: str) -> str:
        bad_chars = ['-', '[', ']']
        for bad_char in bad_chars:
            name = name.replace(bad_char, '')

        return name

    return 'Molecules\t{0}\t{1}'.format(clean(observable.name), _str_from_complex(observable.complex))


def _str_from_rule(rule: Rule, index: int) -> str:
    return '# Rule {5}. rxn: {3}, quant_cont: {4}\n{0} -> {1}   {2}\n'.format(' + '.join(_str_from_complex(x) for x in rule.lhs),
                                                   ' + '.join(_str_from_complex(x) for x in rule.rhs),
                                                   rule.rate.name if rule.rate.name else rule.rate.value,
                                                   str(rule.parent_reaction), str(rule.quant_cont), index + 1)
