from ClickReaction.BaseReaction import BaseReaction, Reactant, Reactants


class AlkalineEsterHydrolysis(BaseReaction):
    """
    Hydrolysis of typical esters that can be removed under alkaline conditions.

    Hydrolyses:
        - methyl esters
        - ethyl esters

    Does not hydrolyse:
        - benzyl esters
        - tBut esters

    amine-boc -> amine
    """

    reactant_names = ["ester"]

    def __init__(self, ester: Reactant):
        self.set_reactants({
            "ester": ester,
        })

    def __runReaction__(self, reactants: Reactants):
        return self._rdReaction.RunReactants((reactants["ester"], ))


smarts = """
    [CX3:1](=[OX1:2])
    -[OX2:3]
    -[$([CH3]),$([CH2]-[CH3])]

    >>

    [C:1](=[O:2])-[O:3]
"""

AlkalineEsterHydrolysis.set_reaction_smarts(smarts)