
import rxncon.semantics.molecule as mol
import rxncon.simulation.SBML.sbml_model as sm

from typing import List, Optional, Set
from libsbml import *


class SBMLBuilder:
    def __init__(self, sbml_storage : sm.SbmlStorage, level = 2, version = 4):
        self.sbml_storage = sbml_storage
        try:
            self.namespace = SBMLNamespaces(level, version)
            self.document = SBMLDocument(self.namespace)
        except ValueError:
            raise SystemExit("SBML Document creation failed")
        self.model = self.document.createModel()

    def process_reaction(self, reactions : List['Reaction']):

        for reaction in reactions:
            sbml_reaction = self.model.createReaction()

            # id is a string "r_" + reactant IDs + "__" +  modifier IDs + "__" + product IDs
            rid="r"
            ["_".join([rid, reactant])for reactant in reaction.reactants]
            rid += "_"
            ["_".join([rid, modifier])for modifier in reaction.modifiers]
            rid += "_"
            ["_".join([rid, product])for product in reaction.products]
            sbml_reaction.setId(rid)

            name = "unkown"
            if reaction.rtype == "1.1.1.1":
                names= "phos"
            elif reaction.rtype == "2.1.1.1":
                names = "ppi"
            sbml_reaction.setName

            sbml_reaction.setReversible(len(reaction.rates) > 1) # if rates holds two rates it has to be reversible
            for reactant in reaction.reactants:
                reactRef = sbml_reaction.createReactant()
                reactRef.setSpecies(reactant)

            for modifier in reaction.modifiers:
                modRef = sbml_reaction.createModifier()
                modRef.setSpecies(modifier)

            for product in reaction.products:
                prodRef =  sbml_reaction.createProduct()
                prodRef.setSpecies(product)

            law = ""
            if reaction.rates[0] is not None:
                law = reaction.rates[0] +" * " + " * ".join(reaction.reactants) + " * " + " * ".join(reaction.modifiers)
                par = self.model.createParameter()
                par.setId(reaction.rates[0])
            if reaction.rates[1] is not None:
                law += " - " + reaction.rates[1] + " * " + " * ".join(reaction.products)
                par = self.model.createParameter()
                par.setId(reaction.rates[1])

            if len(law) > 0:
                kineticLaw = sbml_reaction.createKineticLaw()
                kineticLaw.setMath(parseL3Formula(law))


    def process_species(self, species : List['Species']):
        for s in species:
            sbml_species = self.model.createSpecies()

            sbml_species.setId(s.id.species_id)

            if len(s.molecules) == 1:
                sbl_species.setName(s.molecules[0].name)
            else:
                name = ""
                for molecule in s.molecules:
                    name += molecule.name

            if s.localisation is not None:
                sbml_species.setCompartment(s.localisation)
            else:
                sbml_species.setCompartment("default")

            if s.initial_amount is not None:
                sbml_species.setInitialAmount(s.initial_amount)
            else:
                sbml_species.setInitialAmount(1)     # default

    def process_compartments(self, compartments: List['Compartment']):
        for compartment in compartments:
            c = self.model.createCompartment()
            c.setId(compartment.compartment_name)
            c.setSize(1)

    def build_model(self):
        self.process_compartments(self.sbml_storage.list_of_compartments)
        self.process_species(self.sbml_storage.list_of_species)