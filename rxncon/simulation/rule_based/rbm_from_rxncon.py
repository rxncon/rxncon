import rxncon.core.rxncon_system as rxs
import rxncon.semantics.molecule_from_rxncon as mfr
import rxncon.venntastic.sets as venn


class RuleBasedModelSupervisor:
    def __init__(self, rxncon: rxs.RxnConSystem):
        self.rxncon = rxncon
        self.mol_def_supervisor = mfr.MoleculeDefinitionSupervisor(self.rxncon)

        self.molecules = self.mol_def_supervisor.molecules
        self.rules = []

        self._generate_rules()

    def _generate_rules(self):
        for reaction in self.rxncon.reactions:
            strict_contingency_state_set = mfr.set_of_states_from_contingencies(self.rxncon.strict_contingencies_for_reaction(reaction))
            source_state_set = mfr.source_set_of_states_from_reaction(reaction)

            for molecule in self.molecules:
                instance_set_from_contingencies = mfr.set_of_instances_from_molecule_def_and_set_of_states(
                    self.mol_def_supervisor.molecule_definition_for_name(molecule),
                    strict_contingency_state_set
                )

                disjunct_instance_sets_from_contingencies = venn.gram_schmidt_disjunctify(
                    instance_set_from_contingencies.to_union_list_form()
                )

                instance_set_from_source = mfr.set_of_instances_from_molecule_def_and_set_of_states(
                    self.mol_def_supervisor.molecule_definition_for_name(molecule),
                    source_state_set
                )

                for disjunct_set in disjunct_instance_sets_from_contingencies:
                    lhs = venn.Intersection(venn.Complement(instance_set_from_source), disjunct_set)
                    rhs = venn.Intersection(instance_set_from_source, disjunct_set)




