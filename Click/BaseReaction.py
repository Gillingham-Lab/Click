from typing import Dict, List, Iterable
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem import AllChem

from . import Exceptions


Reactant = Mol
Reactants = Dict[str, Reactant]


class BaseReaction:
    """
    Abstract reaction class to provide common implementations for all reactions.
    """

    _reactants: Reactants
    _smarts: str
    _rdReaction: ChemicalReaction

    def __runReaction__(self, reactants: Reactants) -> List[List[Mol]]:
        """
        Returns all products of all product sets.

        :return: A list of all product sets.
        """
        raise NotImplementedError("You must implement __runReaction__")

    @classmethod
    def set_reaction_smarts(cls, smarts: str) -> None:
        """
        Sets the reaction smarts and creates the rdkit reaction from it.

        :param smarts: A smarts string. All whitespace will be removed.
        """
        cls._smarts = re.sub(r'\s+', '', smarts)
        cls._rdReaction = AllChem.ReactionFromSmarts(cls._smarts)

    def set_reactants(self, reactants: Reactants):
        """
        Sets the reactants.

        :param reactants: A dictionary where each key-value pair associates a reactant name with the corresponding mol.
        :return:

        Example:
            AmideCoupling.set_reactants({"amine": amine_molecule, "acid": acid_molecule})
        """
        self._reactants = reactants

    def get_reactants(self) -> Reactants:
        """
        Returns the reactants as a dictionary, where each key-value pair associates a reactant name with the
        corresponding mol. See set_reactants.

        :return:
        """
        return self._reactants

    def get_products(self, symmetrical_as_one: bool = False) -> List[Mol]:
        """
        Returns a list of all possible products.

        :param symmetrical_as_one: Set to true if symmetrical products should get reduced to one.
        :return:
        """
        productSets = self.__runReaction__(self.get_reactants())

        if productSets is None:
            raise Exception("No product set was returned.")
        elif len(productSets) == 0:
            raise Exceptions.NoProductError("Reaction {} gave no product.".format(type(self)))

        # Retrieve first product of all product sets
        products = []
        productSmiles = []
        for p in productSets:
            # Sanitize product from reaction fragments.
            AllChem.SanitizeMol(p[0])

            if symmetrical_as_one:
                smiles = AllChem.MolToSmiles(p[0])

                if smiles in productSmiles:
                    continue
                else:
                    productSmiles.append(smiles)

            products.append(p[0])

        return products

    def get_product(self, symmetrical_as_one: bool = False) -> Mol:
        """
        Returns one product and raises an exception if multiple products are possible.

        :param symmetrical_as_one: Set to true to remove all but one instance of identical products.
        :return:
        """

        # Get all possible products
        products = self.get_products(symmetrical_as_one=symmetrical_as_one)

        # More than one product is unexpected, raise an error to make the user aware.
        if len(products) > 1:
            raise Exceptions.AmbiguousProductError("Reaction {} gave more than one product sets: {}".format(type(self), [Chem.MolToSmiles(x) for x in products]))

        return products[0]