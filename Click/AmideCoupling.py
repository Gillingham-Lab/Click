import re
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from .Reaction import Reaction, Reactants


class AmideCoupling(Reaction):
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

    def __init__(self, amine: Mol, acid: Mol):
        self.setReactants({
            "amine": amine,
            "acid": acid,
        })

    def __runReaction__(self, reactants: Reactants):
        return self._rdReaction.RunReactants((reactants["amine"], reactants["acid"]))



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