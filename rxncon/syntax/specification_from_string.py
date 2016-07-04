import re
from collections import OrderedDict
import rxncon.core.specification as com
from enum import unique
from rxncon.util.utils import OrderedEnum
import typecheck as tc


@unique
class SpecificationSuffix(OrderedEnum):
    mrna = "mRNA"
    gene = "gene"
    protein = ""


mapping_suffix_to_specification = OrderedDict([(SpecificationSuffix.mrna, com.RnaSpecification),
                                               (SpecificationSuffix.gene, com.DnaSpecification),
                                               (SpecificationSuffix.protein, com.ProteinSpecification),])


def create_specification_from_name_suffix(name, structured_index, domain, subdomain, residue):
    for suffix in mapping_suffix_to_specification:
        if name.endswith(suffix.value):
            name = name[:len(name)-len(suffix.value)]
            return mapping_suffix_to_specification[suffix](name, structured_index, com.DomainResolution(domain, subdomain, residue))


def domain_resolution_from_string(full_domain_string):
    DOMAIN_SUBDOMAIN_RESIDUE_REGEX = '^[\w:-]+\/[\w:-]+\([\w:-]+\)$'
    DOMAIN_RESIDUE_REGEX = '^[\w:-]+\([\w:-]+\)$'
    DOMAIN_SUBDOMAIN_REGEX = '^[\w:-]+\/[\w:-]+$'
    RESIDUE_REGEX = '^\([\w:-]+\)$'
    DOMAIN_REGEX = '^[\w:-]+$'

    if re.match(DOMAIN_SUBDOMAIN_RESIDUE_REGEX, full_domain_string):
        domain = full_domain_string.split('/')[0]
        subdomain = full_domain_string.split('/')[1].split('(')[0]
        residue = full_domain_string.split('/')[1].split('(')[1].strip(')')

    elif re.match(DOMAIN_RESIDUE_REGEX, full_domain_string):
        domain = full_domain_string.split('(')[0]
        subdomain = None
        residue = full_domain_string.split('(')[1].strip(')')

    elif re.match(DOMAIN_SUBDOMAIN_REGEX, full_domain_string):
        domain = full_domain_string.split('/')[0]
        subdomain = full_domain_string.split('/')[1]
        residue = None

    elif re.match(RESIDUE_REGEX, full_domain_string):
        domain = None
        subdomain = None
        residue = full_domain_string.strip('()')

    elif re.match(DOMAIN_REGEX, full_domain_string):
        domain = full_domain_string
        subdomain = None
        residue = None

    else:
        raise SyntaxError('Could not parse specification string {}'.format(specification_string))

    return domain, subdomain, residue


@tc.typecheck
def specification_from_string(specification_string: str) -> com.Specification:
    DOMAIN_DELIMITER = '_'
    SPECIFICATION_DELEMITER = '.'
    #specification_items = specification_string.split('_')
    #assert len(specification_items) <= 2

    items = specification_string.split(DOMAIN_DELIMITER, maxsplit=1)
    name, structured_index = items[0].split('@')

    if len(items) == 1:
        return create_specification_from_name_suffix(name, structured_index, None, None, None)

    elif len(items) == 2:
        full_domain_string = items[1].strip('[]')
        domain, subdomain, residue = domain_resolution_from_string(full_domain_string)
        return create_specification_from_name_suffix(name, structured_index, domain, subdomain, residue)

    else:
        raise SyntaxError('Could not parse specification string {}'.format(specification_string))

