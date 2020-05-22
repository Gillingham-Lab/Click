from Click.BaseReaction import BaseReaction, Reactant, Reactants


class SulfonAmideFormation(BaseReaction):
    """
    Formation of sulfon amides from sulfonyl halogenides and amines

    amine + sulfonylhalogenide -> sulfonamide
    """

    def __init__(self, amine: Reactant, sulfonylhalogenide: Reactant):
        self.set_reactants({
            "amine": amine,
            "sulfonylhalogenide": sulfonylhalogenide,
        })

    def __runReaction__(self, reactants: Reactants):
        return self._rdReaction.RunReactants((reactants["amine"], reactants["sulfonylhalogenide"]))


smarts = """
    [
        $([NX3H3]),
        $([NX4H4]),
        $([NX3H2]-[#6]),
        $([NX4H3]-[#6]),
        $([#7X3H1]([#6])([#6])),
        $([#7X4H2]([#6])([#6]))
    :1]

    .

    [#6:2]-[$([SX4](~[OX1H0])(~[OX1H0]))]-[F,Cl,Br,I;X1]

    >>

    [#7+0:1]-[S+2](-[O-])(-[O-])-[*:2]
"""

SulfonAmideFormation.set_reaction_smarts(smarts)