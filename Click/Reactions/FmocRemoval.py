from Click.BaseReaction import BaseReaction, Reactant, Reactants


class FmocRemoval(BaseReaction):
    """
    Removal of Fmoc from amines

    amine-fmoc -> amine
    """

    def __init__(self, fmocamine: Reactant):
        self.set_reactants({
            "fmocamine": fmocamine,
        })

    def __runReaction__(self, reactants: Reactants):
        return self._rdReaction.RunReactants((reactants["fmocamine"], ))


smarts = """
    [#7:1]
    -C(=O)
    -O
    -C
    -[$([C]1[cR2][cR2][cR2][cR2]1)]

    >>

    [*:1]
"""

FmocRemoval.set_reaction_smarts(smarts)