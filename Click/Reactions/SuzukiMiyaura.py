from Click.BaseReaction import BaseReaction, Reactant, Reactants


class SuzukiMiyaura(BaseReaction):
    """
    Formation of C-C bonds from a boronic acid or ester and an aromatic halogenide (except fluorine).

    boronic acid + arylhalogenide -> C-C bond formation
    boronic ester + arylhalogenide -> C-C bond formation
    trifluoroborates + arylhalogenide -> C-C bond formation
    """

    def __init__(self, boronate: Reactant, halogenide: Reactant):
        self.set_reactants({
            "boronate": boronate,
            "halogenide": halogenide,
        })

    def __runReaction__(self, reactants: Reactants):
        return self._rdReaction.RunReactants((reactants["boronate"], reactants["halogenide"]))


smarts = """
    [#6:1]-[$([B](-O)(-O)),$([B](-F)(-F)(-F))]
    .

    [c:2]-[I,Br,Cl]

    >>

    [*:1]-[*:2]
"""

SuzukiMiyaura.set_reaction_smarts(smarts)