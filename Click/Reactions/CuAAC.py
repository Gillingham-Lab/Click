import re
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from Click.BaseReaction import BaseReaction, Reactants


class CuAAC(BaseReaction):
    """
    Copper-catalyzed alkyne azide cycloaddition to give the 1,4 regioisomer.

    Alkyne + Azide -> 1,4-triazole

    Attributes
    ----------
    alkyne: rdkit.Chem.rdchem.Mol
        A terminal alkyne, or an alkyne iodide.
    azide: rdkit.Chem.rdchem.Mol
        An azide
    """

    def __init__(self, alkyne: Mol, azide: Mol):
        self.setReactants({
            "alkyne": alkyne,
            "azide": azide,
        })

    def __runReaction__(self, reactants: Reactants):
        return self._rdReaction.RunReactants((reactants["alkyne"], reactants["azide"]))


smarts = """
    [C:1]#[$([CH1]),$(C-[I]):2]

    .

    [$([#6]-N=[N+]=[-N]),$([#6]-[N-]-[N+]#N):3]-N~N~N

    >>

    [*:3]n1[c:2][c:1]nn1
"""

CuAAC.setReactionSmarts(smarts)