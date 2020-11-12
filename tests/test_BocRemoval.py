from ClickReaction.Reactions import BocRemoval
from tests import TestHelper


class BocRemovalTest(TestHelper.ReactionTestCase):
    reactant_names = ["bocamine"]

    def setUp(self):
        self.set_reaction(BocRemoval)

    def test_one_product(self):
        reactants = [
            # Boc from primary aliphatic amine
            (("CNC(OC(C)(C)C)=O", ), "CN"),
            # Boc from primary aniline
            (("CC(OC(NC1=CC=CC=C1)=O)(C)C", ), "NC1=CC=CC=C1"),
            # Boc from secondary aliphatic amines
            (("CN(C)C(OC(C)(C)C)=O", ), "CNC"),
            # Boc from secondary anilines
            (("CC(OC(N(C1=CC=CC=C1)C2=CC=CC=C2)=O)(C)C", ), "C1(NC2=CC=CC=C2)=CC=CC=C1"),
            # Boc from mixed anilines / aliphatic amines
            (("CN(C(OC(C)(C)C)=O)C1=CC=CC=C1", ), "CNC1=CC=CC=C1"),
        ]

        self._test_one_product(reactants, self.reactant_names)

    # Those tests should not give any product.
    def test_no_product(self):
        reactants = [
            # Fmoc
            ("CNC(OCC1C(C=CC=C2)=C2C3=C1C=CC=C3)=O", ),

            # Benzyl
            ("CNCC1=CC=CC=C1", ),

            # Dimethylamine
            ("CNC",),

            # O-Methylcarbamates
            ("CN(C(OC)=O)C", ),

            # Ureas
            ("CN(C(NC)=O)C", ),
        ]

        self._test_no_product(reactants, self.reactant_names)
