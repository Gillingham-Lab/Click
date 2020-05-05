import re
from typing import List
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from . import Exceptions


class AmideCoupling:
    """
    Amide coupling reaction: Amine + Carboxylic acid -> Amide

    Attributes
    ----------
    amine: rdkit.Chem.rdchem.Mol
        The amine
    acid: rdkit.Chem.rdchem.Mol
        The carboxylic acid
    """
    _smarts = ""
    _rdReaction = None

    reactants = {}
    product = None

    def __init__(self, amine, acid):
        self.reactants = {
            "amine": amine,
            "acid": acid,
        }

    def __runReaction__(self, amine, acid):
        return self._rdReaction.RunReactants((amine, acid))

    def getProducts(self, symmetrical_as_one: bool = False) -> List[Mol]:
        """
        Returns a list of all possible products.

        Parameters
        ----------
        symmetrical_as_one: bool
            Set to true to remove all but one instance of identical products.

        """
        productSets = self.__runReaction__(**self.reactants)

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

#[
#    $([NH3]),
##    $([NX3H2]-[CX4]),
#    $([NX3H1](-[CX4])(-[CX4]))
#:1]

AmideCoupling._smarts = re.sub(r'\s+', '', """
    [
        $([NX3H3]),
        $([NX4H4]),
        $([NX3H2]-[CX4]),
        $([NX4H3]-[CX4]),
        $([NX3H1](-[CX4])(-[CX4])),
        $([NX4H2](-[CX4])(-[CX4]))
    :1]

    .
    
    [C:2](=[OX1:3])-[$([OX2H1]),$([O-X1])]
    
    >>
    
    [*+0:1]-[*:2](=[*:3])
""")
AmideCoupling._rdReaction = AllChem.ReactionFromSmarts(AmideCoupling._smarts)