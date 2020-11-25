from ClickReaction.Reactions.AmideCoupling import AmideCoupling
from ClickReaction.Reactions.AmideCouplingWithAnilines import AmideCouplingWithAnilines
from ClickReaction.Reactions.CuAAC import CuAAC
from ClickReaction.Reactions.SulfonAmideFormation import SulfonAmideFormation

# Cross couplings
from ClickReaction.Reactions.SuzukiMiyaura import SuzukiMiyaura

# Protecting group removals
from ClickReaction.Reactions.FmocRemoval import FmocRemoval
from ClickReaction.Reactions.BocRemoval import BocRemoval
from ClickReaction.Reactions.AlkalineEsterHydrolysis import AlkalineEsterHydrolysis

# All reactions as a dictionary.
all_reactions = {
    "AmideCoupling": AmideCoupling,
    "AmideCouplingWithAnilines": AmideCouplingWithAnilines,
    "CuAAC": CuAAC,
    "SulfonAmideFormation": SulfonAmideFormation,
    "SuzukiMiyaura": SuzukiMiyaura,
    "FmocRemoval": FmocRemoval,
    "BocRemoval": BocRemoval,
    "AlkalineEsterHydrolysis": AlkalineEsterHydrolysis,
}
