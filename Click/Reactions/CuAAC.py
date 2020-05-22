from Click.BaseReaction import BaseReaction, Reactant, Reactants


class CuAAC(BaseReaction):
    """
    Copper-catalyzed alkyne azide cycloaddition to give the 1,4 regioisomer.

    alkyne + azide -> 1,4-triazole
    """

    def __init__(self, alkyne: Reactant, azide: Reactant):
        self.set_reactants({
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

CuAAC.set_reaction_smarts(smarts)