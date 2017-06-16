"""Module containing the functions boolnet_from_boolean_model, boolnet_strs_from_rxncon and
the class QuantitativeContingencyStrategy."""

from enum import Enum
from typing import Tuple, Dict

from rxncon.core.rxncon_system import RxnConSystem
from rxncon.simulation.boolean.boolean_model import BooleanModel, Target, ReactionTarget, KnockoutTarget, \
    OverexpressionTarget, StateTarget, UpdateRule, SmoothingStrategy, KnockoutStrategy, OverexpressionStrategy, \
    boolean_model_from_rxncon
from rxncon.venntastic.sets import Set as VennSet, ValueSet, Complement, Intersection, Union, EmptySet, UniversalSet


def boolnet_from_boolean_model(boolean_model: BooleanModel) -> Tuple[str, Dict[str, str], Dict[str, bool]]:
    """Translates the boolean model into BoolNet syntax.

    Returns:
        1. The boolean model in BoolNet syntax,
        2. The (BoolNet name, rxncon name) mapping, since BoolNet is picky about characters,
        3. The initial conditions."""

    def initialize_boolnet_names() -> None:
        nonlocal boolnet_names
        nonlocal reaction_index
        nonlocal state_index
        nonlocal knockout_index
        nonlocal overexpression_index
        for update_rule in boolean_model.update_rules:
            target = update_rule.target
            if target in boolnet_names.keys():
                raise AssertionError
            elif isinstance(target, ReactionTarget):
                name = 'R{}'.format(reaction_index)
                boolnet_names[target] = name
                reaction_index += 1
            elif isinstance(target, KnockoutTarget):
                name = 'K{}'.format(knockout_index)
                boolnet_names[target] = name
                knockout_index += 1
            elif isinstance(target, OverexpressionTarget):
                name = 'O{}'.format(overexpression_index)
                boolnet_names[target] = name
                overexpression_index += 1
            elif isinstance(target, StateTarget):
                name = 'S{}'.format(state_index)
                boolnet_names[target] = name
                state_index += 1
            else:
                raise AssertionError

    def str_from_factor(factor: VennSet) -> str:
        if isinstance(factor, ValueSet):
            return boolnet_names[factor.value]
        elif isinstance(factor, Complement):
            return '!({})'.format(str_from_factor(factor.expr))
        elif isinstance(factor, Intersection):
            return '({})'.format(' & '.join(str_from_factor(x) for x in factor.exprs))
        elif isinstance(factor, Union):
            return '({})'.format(' | '.join(str_from_factor(x) for x in factor.exprs))
        elif isinstance(factor, EmptySet):
            return '0'
        elif isinstance(factor, UniversalSet):
            return '1'
        else:
            raise AssertionError('Could not parse factor {}'.format(factor))

    def str_from_update_rule(update_rule: UpdateRule) -> str:
        return '{0}, {1}'.format(boolnet_names[update_rule.target],
                                 str_from_factor(update_rule.factor))

    boolnet_names = {}  # type: Dict[Target, str]

    reaction_index = 0
    state_index = 0
    knockout_index = 0
    overexpression_index = 0

    initialize_boolnet_names()

    def sort_key(rule_str: str) -> Tuple[str, int]:
        target = rule_str.split(',')[0].strip()
        return target[0], int(target[1:])

    rule_strs = sorted([str_from_update_rule(x) for x in boolean_model.update_rules], key=sort_key)

    return 'targets, factors\n' + '\n'.join(rule for rule in rule_strs) + '\n', \
           {name: str(target) for target, name in boolnet_names.items()}, \
           {boolnet_names[target]: value for target, value in boolean_model.initial_conditions.target_to_value.items()}


class QuantitativeContingencyStrategy(Enum):
    """Strategy for dealing with quantitative contingencies k+ and k-, ignore them or make them into ! resp. x."""
    strict = 'strict'
    ignore = 'ignore'


def boolnet_strs_from_rxncon(rxncon: RxnConSystem, smoothing_strategy: SmoothingStrategy,
                             knockout_strategy: KnockoutStrategy,
                             overexpression_strategy: OverexpressionStrategy,
                             k_plus_strategy: QuantitativeContingencyStrategy,
                             k_minus_strategy: QuantitativeContingencyStrategy) \
        -> Tuple[str, str, str]:
    """Returns a triple of strs:
         1. The BoolNet model,
         2. The mapping between BoolNet names and rxncon names, sorted by BoolNet name, and
         3. The initial values, sorted by BoolNet name.
    Different from boolnet_from_boolean_model: first converts from rxncon model (and therefore needs these
    strategies) and also returns strings that can directly be written to a file."""
    def sort_key(key_val_pair):
        k, v = key_val_pair
        return k[0], int(k[1:])

    if k_plus_strategy == QuantitativeContingencyStrategy.strict:
        k_plus_strict = True
    elif k_plus_strategy == QuantitativeContingencyStrategy.ignore:
        k_plus_strict = False
    else:
        raise AssertionError('Unknown QuantitativeContingencyStrategy {}'.format(k_plus_strategy))

    if k_minus_strategy == QuantitativeContingencyStrategy.strict:
        k_minus_strict = True
    elif k_minus_strategy == QuantitativeContingencyStrategy.ignore:
        k_minus_strict = False
    else:
        raise AssertionError('Unknown QuantitativeContingencyStrategy {}'.format(k_minus_strategy))

    model_str, symbol_dict, initial_val_dict = \
        boolnet_from_boolean_model(boolean_model_from_rxncon(rxncon, smoothing_strategy=smoothing_strategy,
                                                             knockout_strategy=knockout_strategy,
                                                             overexpression_strategy=overexpression_strategy,
                                                             k_plus_strict=k_plus_strict,
                                                             k_minus_strict=k_minus_strict))

    symbol_str = '\n'.join('{0}, {1}'.format(boolnet_sym, rxncon_sym) for boolnet_sym, rxncon_sym
                           in sorted(symbol_dict.items(), key=sort_key)) + '\n'

    initial_val_str = '\n'.join('{0}, {1: <5}  , #  {2}'.format(boolnet_sym, initial_val, symbol_dict[boolnet_sym])
                                for boolnet_sym, initial_val in sorted(initial_val_dict.items(), key=sort_key)) + '\n'

    return model_str, symbol_str, initial_val_str
