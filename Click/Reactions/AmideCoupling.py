from Click.BaseReaction import BaseReaction, Reactant, Reactants


class AmideCoupling(BaseReaction):
    """
    Amide coupling reaction to form amides.

    amine + carboxylic acid -> Amide
    """

    def __init__(self, amine: Reactant, acid: Reactant):
        self.set_reactants({
            "amine": amine,
            "acid": acid,
        })

    def __runReaction__(self, reactants: Reactants):
        return self._rdReaction.RunReactants((reactants["amine"], reactants["acid"]))


smarts = """
    [
        $([NX3H3]),
        $([NX4H4]),
        $([NX3H2]-[CX4]),
        $([NX4H3]-[CX4]),
        $([NX3H1](-[CX4])(-[CX4])),
        $([NX4H2](-[CX4])(-[CX4]))
    :1]

    .
    
    [C:2]
        (=[OX1:3])
        -[
            $([OX2]-N1C(=O)CCC1(=O)),
            $([OX2]-c1c(-F)c(-F)c(-F)c(-F)c1(-F)),
            $([OX2]-c1ccc(-N(~O)(~O))cc1),
            $([OX2H1]),
            $([O-X1])
        ]
    
    >>
    
    [*+0:1]-[*:2](=[*:3])
"""

AmideCoupling.set_reaction_smarts(smarts)