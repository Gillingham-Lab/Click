from typing import Dict, List, Set
import re
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

    def __runReaction__(self, reactants: Reactants) -> Set[Set[Mol]]:
        """
        Returns all products of all product sets.

        :return:
        """
        raise NotImplementedError("You must implement __runReaction__")

    @classmethod
    def setReactionSmarts(cls, smarts: str) -> None:
        """
        Sets a reaction smarts and creates the rdkit reaction.

        :param smarts: A smarts string. All whitespace will be removed.
        :return:
        """
        cls._smarts = re.sub(r'\s+', '', smarts)
        cls._rdReaction = AllChem.ReactionFromSmarts(cls._smarts)

    def setReactants(self, reactants: Reactants):
        """
        Sets the reactants.

        :param reactants:
        :return:
        """
        self._reactants = reactants

    def getReactants(self) -> Reactants:
        """
        Returns the reactants as a dictionary.

        :return:
        """
        return self._reactants

    def getProducts(self, symmetrical_as_one: bool = False) -> List[Mol]:
        """
        Returns a list of all possible products.

        Parameters
        ----------
        symmetrical_as_one: bool
            Set to true to remove all but one instance of identical products.

        """
        productSets = self.__runReaction__(self.getReactants())

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

    def getProduct(self, symmetrical_as_one: bool = False) -> Mol:
        """
        Returns one product and raises an exception if multiple products are possible.

        Parameters
        ----------
        symmetrical_as_one: bool
            Set to true to remove all but one instance of identical products.
            This allows reactions that create two or more identical products because of symmetrical a symmetric reactant
             to only return one product.

        """

        # Get all possible products
        products = self.getProducts(symmetrical_as_one=symmetrical_as_one)

        # More than one product is unexpected, raise an error to make the user aware.
        if len(products) > 1:
            raise Exceptions.AmbiguousProductError("Reaction {} gave more than one product sets.".format(type(self)))

        return products[0]