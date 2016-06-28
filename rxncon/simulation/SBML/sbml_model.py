from typing import List, Optional, Set
import rxncon.semantics.molecule as mol


class SbmlStorage:
    def __init__(self, species: List['Species'], reactions: List['Reaction'], compartments : List['Compartment']):
        self.list_of_species = species
        self.list_of_reactions = reactions
        self.list_of_compartments = compartments


class Species:
    pass


class SingleSpecies(Species):
    def __init__(self, molecules: List[mol.Molecule], inital_amount: Optional['InitialAmount'], id: 'SpeciesID'):
        self.list_of_molecules = molecules
        self.id = id
        self.inital_amount = inital_amount
        self._validate()
        # name is different in SBML and CellDesigner through the molecule description level
        # someway to make clear the structure of an complex , maybe subcomplexes or whatever?

    def _validate(self):
        assert len(self.listOfMolecules) == 1


class ComplexSpecies(Species):
    def __init__(self, molecules: List[mol.Molecule], inital_amount: Optional['InitialAmount'], id: 'SpeciesID'):
        self.list_of_molecules = molecules
        self.id = id
        self.inital_amount = inital_amount
        self._validate()
    def _validate(self):
        assert len(self.listOfMolecules) > 1


class Reaction:
    def __init__(self, reactantRef : Set[SpeciesID], modRef : Set[SpeciesID], prodRef : Set[SpeciesID], type : str, rates : [int]):
        self.reactants = reactantRef    # reference to the list of molecules
        self.modifiers = modRef
        self.products = prodRef
        self.rtype = type
        self.rates = rates  # list of rates, either one or two, if two the reaction is reversible


class SpeciesID:
    def __init__(self, species_id: str):
        self.species_id = species_id
        self._validate()
    def _validate(self):
        pass


class Compartment:
    def __init__(self, compartment_name: str):
        self.compartment_name = compartment_name


class Compartments:
    def __init__(self, valid_compartments: Set['Compartment']):
        self.compartments = valid_compartments


class InitialAmount:
    def __init__(self, amount: int):
        self.amount = amount


class SpeciesReference:
    pass


class Reactant(SpeciesReference):
    pass


class Product(SpeciesReference):
    pass


class Modifier(SpeciesReference):
    pass