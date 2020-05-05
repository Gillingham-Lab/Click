from rdkit.Chem import AllChem

from . import Exceptions


class AmideCoupling():
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

    def getProduct(self):
        productSets =  self.__runReaction__(**self.reactants)

        if productSets is None:
            raise Exception("No product set was returned.")
        elif len(productSets) == 0:
            raise Exceptions.NoProductError("Reaction {} gave no product.".format(type(self)))
        elif len(productSets) > 1:
            raise Exceptions.AmbiguousProductError("Reaction {} gave more than one product sets.".format(type(self)))

        product = productSets[0][0]

        AllChem.SanitizeMol(product)
        return product


AmideCoupling._smarts = "[$([NH3]),$([NX3H2]-[CX4]),$([NX3H1](-[CX4])(-[CX4])):1].[C:2](=[OX1:3])-[OH1]>>[N:1]-[C:2](=[O:3])"
AmideCoupling._rdReaction = AllChem.ReactionFromSmarts(AmideCoupling._smarts)