from rxncon.core.rxncon_system import RxnConSystem
from rxncon.semantics.molecule import MoleculeDef


def molecule_defs_from_rxncon(rxnconsys: RxnConSystem):
    class MolDefDict(dict):
        def __getitem__(self, key):
            if key not in self.keys():
                self[key] = MoleculeDef(key)

            return dict.__getitem__(self, key)

    molecule_defs = MolDefDict()

    for state in rxnconsys.product_states + rxnconsys.source_states:
        for spec in state.specs:
            molecule_defs[spec.to_component_spec()].add_state(state)

    return molecule_defs



