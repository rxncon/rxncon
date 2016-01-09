from collections import namedtuple
from enum import Enum, unique
from typing import List

import rxncon.core.component as com
import rxncon.core.error as err
import rxncon.core.state as sta


@unique
class Verb(Enum):
    phosphorylation             = 'P+'
    dephosphorylation           = 'P-'
    autophosphorylation         = 'AP'
    phosphotransfer             = 'PT'
    guanine_nucleotide_exchange = 'GEF'
    gtpase_activation           = 'GAP'
    ubiquination                = 'Ub+'
    proteolytic_cleavage        = 'CUT'
    protein_protein_interaction = 'ppi'
    intra_protein_interaction   = 'ipi'
    non_protein_interaction     = 'i'
    binding_to_dna              = 'BIND'
    degradation                 = 'DEG'
    synthesis                   = 'SYN'


@unique
class ReactionClass(Enum):
    covalent_modification = 1
    interaction           = 2
    synthesis_degradation = 3
    translocation         = 4


@unique
class Directionality(Enum):
    reversible   = 1
    irreversible = 2


# @todo Bidirectional needs separate identifier?
class Influence(Enum):
    positive      = 1
    negative      = 2
    transfer      = 3
    bidirectional = 1


@unique
class Isomerism(Enum):
    undefined = None
    trans     = 1
    cis       = 2


@unique
class CovalentReactionModifier(Enum):
    undefined = None
    phosphor  = 'P'
    ubiquitin = 'Ub'
    truncated = 'Truncated'


VERB_REACTION_TABLE = {
    Verb.phosphorylation:              [ReactionClass.covalent_modification,
                                        Directionality.reversible,
                                        Influence.positive,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.phosphor],
    Verb.dephosphorylation:            [ReactionClass.covalent_modification,
                                        Directionality.reversible,
                                        Influence.negative,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.phosphor],
    Verb.autophosphorylation:          [ReactionClass.covalent_modification,
                                        Directionality.reversible,
                                        Influence.positive,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.phosphor],
    Verb.phosphotransfer:              [ReactionClass.covalent_modification,
                                        Directionality.reversible,
                                        Influence.transfer,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.phosphor],
    Verb.guanine_nucleotide_exchange:  [ReactionClass.covalent_modification,
                                        Directionality.reversible,
                                        Influence.positive,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.phosphor],
    Verb.gtpase_activation:            [ReactionClass.covalent_modification,
                                        Directionality.reversible,
                                        Influence.negative,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.phosphor],
    Verb.ubiquination:                 [ReactionClass.covalent_modification,
                                        Directionality.reversible,
                                        Influence.positive,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.ubiquitin],
    Verb.proteolytic_cleavage:         [ReactionClass.covalent_modification,
                                        Directionality.irreversible,
                                        Influence.positive,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.truncated],
    Verb.protein_protein_interaction:  [ReactionClass.interaction,
                                        Directionality.reversible,
                                        Influence.bidirectional,
                                        Isomerism.trans,
                                        CovalentReactionModifier.undefined],
    Verb.intra_protein_interaction:    [ReactionClass.interaction,
                                        Directionality.reversible,
                                        Influence.bidirectional,
                                        Isomerism.cis,
                                        CovalentReactionModifier.undefined],
    Verb.non_protein_interaction:      [ReactionClass.interaction,
                                        Directionality.reversible,
                                        Influence.positive,
                                        Isomerism.trans,
                                        CovalentReactionModifier.undefined],
    Verb.binding_to_dna:               [ReactionClass.interaction,
                                        Directionality.reversible,
                                        Influence.positive,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.undefined],
    Verb.synthesis:                    [ReactionClass.synthesis_degradation,
                                        Directionality.irreversible,
                                        Influence.positive,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.undefined],
    Verb.degradation:                  [ReactionClass.synthesis_degradation,
                                        Directionality.irreversible,
                                        Influence.negative,
                                        Isomerism.undefined,
                                        CovalentReactionModifier.undefined]
}


class Reaction:
    def __init__(self, full_name: str, subject: com.Component, verb: Verb, object: com.Component,
                 reaction_class: ReactionClass, directionality: Directionality, influence: Influence,
                 isomerism: Isomerism, modifier: CovalentReactionModifier):
        self.full_name = full_name
        self.subject = subject
        self.verb = verb
        self.object = object

        self.reaction_class = reaction_class
        self.directionality = directionality
        self.influence = influence
        self.isomerism = isomerism
        self.modifier = modifier

        self.source = None   # type: Optional[sta.State]
        self.product = None  # type: Optional[sta.State]
        self.classification_code = None  # type: str

        self._determine_source_product_states()
        self._determine_classification_code()

    def __repr__(self) -> str:
        return self.full_name

    def __eq__(self, other: 'Reaction') -> bool:
        assert isinstance(other, Reaction)
        return self.full_name == other.full_name and self.subject == other.subject and self.verb == other.verb and \
            self.object == other.object and self.reaction_class == other.reaction_class and self.directionality == other.directionality and \
            self.influence == other.influence and self.isomerism == other.isomerism and self.modifier == other.modifier

    @property
    def constant_components(self) -> List[com.Component]:
        # @todo implement this?
        return []

    def _determine_source_product_states(self):
        self.source, self.product = states_from_reaction(self)

    def _determine_classification_code(self):
        properties = [self.reaction_class, self.directionality, self.influence]
        if self.isomerism.value:
            properties.append(self.isomerism)

        self.classification_code = '.'.join([str(p.value) for p in properties])


class InputReaction(Reaction):
    def __init__(self, full_name: str):
        self.full_name = full_name

    def __eq__(self, other: Reaction) -> bool:
        assert isinstance(other, Reaction)

        if isinstance(other, InputReaction):
            return self.full_name == other.full_name

        else:
            return False


SourceStateProductState = namedtuple('SourceStateProductState', ['source_state', 'product_state'])


def states_from_reaction(reaction: Reaction) -> SourceStateProductState:
    if isinstance(reaction, InputReaction):
        return SourceStateProductState(None, None)

    elif reaction.reaction_class == ReactionClass.covalent_modification:
        return _covalent_modification_states_from_reaction(reaction)

    elif reaction.reaction_class == ReactionClass.interaction:
        return _interaction_states_from_reaction(reaction)

    elif reaction.reaction_class == ReactionClass.synthesis_degradation:
        return _synthesis_degradation_states_from_reaction(reaction)

    elif reaction.reaction_class == ReactionClass.translocation:
        return _translocation_states_from_reaction(reaction)

    else:
        raise err.RxnConLogicError('Non-exhaustive switch statement in state_from_reaction for case {}'.format(reaction))


def _covalent_modification_states_from_reaction(reaction: Reaction) -> SourceStateProductState:
    if reaction.modifier == CovalentReactionModifier.phosphor:
        modifier = sta.StateModifier.phosphor

    elif reaction.modifier == CovalentReactionModifier.ubiquitin:
        modifier = sta.StateModifier.ubiquitin

    elif reaction.modifier == CovalentReactionModifier.truncated:
        modifier = sta.StateModifier.truncated

    else:
        raise err.RxnConLogicError('Could not map rxn modifier {0} to state modifier'.format(reaction.modifier))

    source = product = None

    if reaction.influence == Influence.positive:
        product = sta.CovalentModificationState(reaction.object, modifier)

    elif reaction.influence == Influence.negative:
        source = sta.CovalentModificationState(reaction.object, modifier)

    elif reaction.influence == Influence.transfer:
        source = sta.CovalentModificationState(reaction.subject, modifier)
        product = sta.CovalentModificationState(reaction.object, modifier)

    else:
        raise err.RxnConLogicError('Could not determine product/source pair for reaction {}'.format(reaction))

    return SourceStateProductState(source, product)


def _interaction_states_from_reaction(reaction: Reaction) -> SourceStateProductState:
    source = None
    product = sta.InteractionState(reaction.subject, reaction.object)

    return SourceStateProductState(source, product)


def _synthesis_degradation_states_from_reaction(reaction: Reaction) -> SourceStateProductState:
    source = product = None

    if reaction.influence == Influence.positive:
        product = sta.SynthesisDegradationState(reaction.object)

    elif reaction.influence == Influence.negative:
        source = sta.SynthesisDegradationState(reaction.object)

    else:
        raise err.RxnConLogicError('Could not determine syn/deg state for reaction {}'.format(reaction))

    return SourceStateProductState(source, product)


def _translocation_states_from_reaction(reaction: Reaction) -> SourceStateProductState:
    source = product = None

    # @todo Map the reaction modifiers to state modifiers.
    assert False

    if reaction.influence == Influence.positive:
        product = sta.TranslocationState(reaction.object, reaction.modifier)

    elif reaction.influence == Influence.negative:
        source = sta.TranslocationState(reaction.object, reaction.modifier)

    else:
        raise err.RxnConLogicError('Could not determine product/source pair for reaction {}'.format(reaction))

    return SourceStateProductState(source, product)