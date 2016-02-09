from collections import namedtuple
from enum import Enum, unique
from typing import List, Optional
import typecheck as tc

import rxncon.core.component as com
import rxncon.core.error as err
import rxncon.core.state as sta
import rxncon.syntax.string_from_rxncon as sfr


@unique
class Verb(Enum):
    phosphorylation             = 'p+'
    dephosphorylation           = 'p-'
    autophosphorylation         = 'ap'
    phosphotransfer             = 'pt'
    guanine_nucleotide_exchange = 'gef'
    gtpase_activation           = 'gap'
    ubiquination                = 'ub+'
    proteolytic_cleavage        = 'cut'
    protein_protein_interaction = 'ppi'
    intra_protein_interaction   = 'ipi'
    non_protein_interaction     = 'i'
    binding_to_dna              = 'bind'
    degradation                 = 'deg'
    synthesis                   = 'syn'


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

# @todo lower case.
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
    @tc.typecheck
    def __init__(self, subject: com.Component, verb: Verb, object: com.Component,
                 reaction_class: ReactionClass, directionality: Directionality, influence: Influence,
                 isomerism: Isomerism, modifier: CovalentReactionModifier):
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

    def __str__(self):
        return sfr.string_from_reaction(self)

    @tc.typecheck
    def __eq__(self, other: 'Reaction') -> bool:
        return self.subject == other.subject and self.verb == other.verb and self.object == other.object and \
            self.reaction_class == other.reaction_class and self.directionality == other.directionality and \
            self.influence == other.influence and self.isomerism == other.isomerism and self.modifier == other.modifier

    @property
    def components(self):
        return [self.subject, self.object]

    def _determine_source_product_states(self):
        self.source, self.product = states_from_reaction(self)

    def _determine_classification_code(self):
        properties = [self.reaction_class, self.directionality, self.influence]
        if self.isomerism.value:
            properties.append(self.isomerism)

        self.classification_code = '.'.join([str(p.value) for p in properties])


class OutputReaction(Reaction):
    def __init__(self, name: str):
        self.name = name

    def __eq__(self, other: Reaction) -> bool:
        assert isinstance(other, Reaction)

        if isinstance(other, OutputReaction):
            return self.name == other.name

        else:
            return False


SourceStateProductState = namedtuple('SourceStateProductState', ['source_state', 'product_state'])


def states_from_reaction(reaction: Reaction) -> SourceStateProductState:
    if isinstance(reaction, OutputReaction):
        return SourceStateProductState(None, None)

    elif reaction.reaction_class == ReactionClass.covalent_modification:
        return _covalent_modification_states_from_reaction(reaction)

    # @todo Is this correct?
    elif reaction.reaction_class == ReactionClass.interaction and \
            (reaction.isomerism == Isomerism.trans or reaction.isomerism == Isomerism.undefined):
        return _inter_protein_interaction_states_from_reaction(reaction)

    elif reaction.reaction_class == ReactionClass.interaction and reaction.isomerism == Isomerism.cis:
        return _intra_protein_interaction_states_from_reaction(reaction)

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


def _inter_protein_interaction_states_from_reaction(reaction: Reaction) -> SourceStateProductState:
    source = None
    product = sta.InterProteinInteractionState(reaction.subject, reaction.object)

    return SourceStateProductState(source, product)


def _intra_protein_interaction_states_from_reaction(reaction: Reaction) -> SourceStateProductState:
    source = None
    product = sta.IntraProteinInteractionState(reaction.subject, reaction.object)

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