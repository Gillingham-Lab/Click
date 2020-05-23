from Click.Reactions import FmocRemoval
from tests import TestHelper


class FmocRemovalTest(TestHelper.ReactionTestCase):
    reactant_names = ["fmocamine"]

    def setUp(self):
        self.set_reaction(FmocRemoval)

    def test_one_product(self):
        reactants = [
            # Fmoc from primary aliphatic amine
            (("CNC(OCC1C(C=CC=C2)=C2C3=C1C=CC=C3)=O", ), "CN"),
            # Fmoc from primary aniline
            (("O=C(NC1=CC=CC=C1)OCC2C(C=CC=C3)=C3C4=C2C=CC=C4", ), "NC1=CC=CC=C1"),
            # Fmoc from secondary aliphatic amines
            (("CN(C)C(OCC1C(C=CC=C2)=C2C3=C1C=CC=C3)=O", ), "CNC"),
            # Fmoc from secondary anilines
            (("O=C(N(C1=CC=CC=C1)C2=CC=CC=C2)OCC3C(C=CC=C4)=C4C5=C3C=CC=C5", ), "C1(NC2=CC=CC=C2)=CC=CC=C1"),
            # Fmoc from mixed anilines / aliphatic amines
            (("CN(C(OCC1C(C=CC=C2)=C2C3=C1C=CC=C3)=O)C4=CC=CC=C4", ), "CNC1=CC=CC=C1"),
        ]

        self._test_one_product(reactants, self.reactant_names)

    # Those tests should not give any product.
    def test_no_product(self):
        reactants = [
            # Boc
            ("CNC(OC(C)(C)C)=O", ),

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
