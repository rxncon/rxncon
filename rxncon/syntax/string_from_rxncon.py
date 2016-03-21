# Note: Here be no type annotations, because this would require importing the state, component and reaction modules,
# which would induce a circular dependency.


def string_from_component(component) -> str:
    if component.domain and component.subdomain and component.residue:
        return '{0}_[{1}/{2}({3})]'.format(component.name, component.domain, component.subdomain, component.residue)

    elif component.domain and not component.subdomain and component.residue:
        return '{0}_[{1}({2})]'.format(component.name, component.domain, component.residue)

    elif component.domain and component.subdomain and not component.residue:
        return '{0}_[{1}/{2}]'.format(component.name, component.domain, component.subdomain)

    elif not component.domain and component.subdomain and component.residue:
        return '{0}_[/{1}({2})]'.format(component.name, component.subdomain, component.residue)

    elif not component.domain and not component.subdomain and component.residue:
        return '{0}_[({1})]'.format(component.name, component.residue)

    elif not component.domain and component.subdomain and not component.residue:
        return '{0}_[/{1}]'.format(component.name, component.subdomain)

    elif component.domain and not component.subdomain and not component.residue:
        return '{0}_[{1}]'.format(component.name, component.domain)

    elif not component.domain and not component.subdomain and not component.residue:
        return '{0}'.format(component.name)

    else:
        raise AssertionError


def string_from_reaction(reaction) -> str:
    return '{0}_{1}_{2}'.format(reaction.subject, reaction.verb.value, reaction.object)


def string_from_inter_protein_interaction_state(state) -> str:
    return '{0}--{1}'.format(state.first_component, state.second_component)


def string_from_intra_protein_interaction_state(state) -> str:
    return '{0}--[{1}]'.format(state.first_component, state.second_component.domain)


def string_from_covalent_modification_state(state) -> str:
    return '{0}-{{{1}}}'.format(state.substrate, state.modifier.value)


def string_from_translocation_state(state) -> str:
    return '{0}-{{{1}}}'.format(state.substrate, state.compartment.value)


def string_from_synthesis_degradation_state(state) -> str:
    return '{}'.format(state.component)
