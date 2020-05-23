from Click.BaseReaction import BaseReaction, Reactant, Reactants


class BocRemoval(BaseReaction):
    """
    Removal of Boc from amines

    amine-boc -> amine
    """

    def __init__(self, bocamine: Reactant):
        self.set_reactants({
            "bocamine": bocamine,
        })

    def __runReaction__(self, reactants: Reactants):
        return self._rdReaction.RunReactants((reactants["bocamine"], ))


smarts = """
    [#7:1]
    -C(=O)
    -O
    -[$(C(-[CH3])(-[CH3])(-[CH3]))]

    >>

    [*:1]
"""

BocRemoval.set_reaction_smarts(smarts)