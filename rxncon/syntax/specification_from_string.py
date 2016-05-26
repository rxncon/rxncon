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


def create_specification_from_name_suffix(name, domain, subdomain, residue, specification_type):
    for suffix in mapping_suffix_to_specification:
        if specification_type == suffix.value:
            return mapping_suffix_to_specification[suffix](name, com.DomainResolution(domain, subdomain, residue))


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
    specification_items = specification_string.split('.')
    assert len(specification_items) <= 2
    if len(specification_items) > 1:
        specification_type = specification_items[1]
    else:
        specification_type = ""
    items = specification_items[0].split(DOMAIN_DELIMITER, maxsplit=1)

    if len(items) == 1:
        return create_specification_from_name_suffix(items[0], None, None, None, specification_type)

    elif len(items) == 2:
        name = items[0]
        full_domain_string = items[1].strip('[]')
        domain, subdomain, residue = domain_resolution_from_string(full_domain_string)
        return create_specification_from_name_suffix(name, domain, subdomain, residue, specification_type)

    else:
        raise SyntaxError('Could not parse specification string {}'.format(specification_string))

