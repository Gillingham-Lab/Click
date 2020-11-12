from ClickReaction.Reactions import SulfonAmideFormation
from tests import TestHelper


class SulfonAmideFormationTest(TestHelper.ReactionTestCase):
    reactant_names = ["amine", "sulfonylhalogenide"]

    def setUp(self):
        self.set_reaction(SulfonAmideFormation)

    def test_one_product(self):
        reactants = [
            # Ammonia and sulfonyl halogenides
            (("N", "F[S+2]([O-])([O-])CC"), "N[S+2]([O-])([O-])CC"),
            (("N", "Cl[S+2]([O-])([O-])CC"), "N[S+2]([O-])([O-])CC"),
            (("N", "Br[S+2]([O-])([O-])CC"), "N[S+2]([O-])([O-])CC"),
            (("N", "I[S+2]([O-])([O-])CC"), "N[S+2]([O-])([O-])CC"),

            # Primary amines and sulfonyl chlorides
            (("CN", "Cl[S+2]([O-])([O-])CC"), "CN[S+2]([O-])([O-])CC"),
            (("C[NH3+]", "Cl[S+2]([O-])([O-])CC"), "CN[S+2]([O-])([O-])CC"),

            # Secondary amine and sulfonyl chloride
            (("CNC", "Cl[S+2]([O-])([O-])CC"), "CN(-C)[S+2]([O-])([O-])CC"),
            (("C[NH2+]C", "Cl[S+2]([O-])([O-])CC"), "CN(-C)[S+2]([O-])([O-])CC"),

            # Primary aniline and sulfonyl chlorides
            (("c1ccccc1N", "Cl[S+2]([O-])([O-])CC"), "c1ccccc1N[S+2]([O-])([O-])CC"),
            (("c1ccccc1[NH3+]", "Cl[S+2]([O-])([O-])CC"), "c1ccccc1N[S+2]([O-])([O-])CC"),

            # Secondary aniline and sulfonyl chloride
            (("c1ccccc1NC", "Cl[S+2]([O-])([O-])CC"), "c1ccccc1N(C)[S+2]([O-])([O-])CC"),

            # Pyrrole or indole
            (("C1=CC=CN1", "Cl[S+2]([O-])([O-])CC"), "c1cccn1[S+2]([O-])([O-])CC"),
            (("C1=CC=CN1", "Cl[S+2]([O-])([O-])CC"), "c1cccn1[S+2]([O-])([O-])CC"),
            (("C12=CC=CC=C1NC=C2", "Cl[S+2]([O-])([O-])CC"), "c12ccccc1n([S+2]([O-])([O-])CC)cc2"),
            (("c1ccc2[nH]ccc2c1", "Cl[S+2]([O-])([O-])CC"), "c12ccccc1n([S+2]([O-])([O-])CC)cc2"),

            # Special test: This is a standard sulfon amide directly from ChemDraw. The tests should be able to digest this.
            (("CN", "ClS(CC)(=O)=O"), "CN[S+2]([O-])([O-])CC"),
        ]

        self._test_one_product(reactants, self.reactant_names)

    # Those tests should not give any product.
    def test_no_product(self):
        reactants = [
            # Ammonium and sulfonic acid
            ("N", "O[S+2]([O-])([O-])CC"),

            # Ammonium and sulfonate
            ("N", "[O-][S+2]([O-])([O-])CC"),

            # No tertiary amine
            ("N(C)(C)(C)", "F[S+2]([O-])([O-])CC"),

            # No tertiary pyrrole or indole
            ("c1cccn(C)1", "F[S+2]([O-])([O-])CC"),
            ("c1ccc2n(C)ccc2c1", "F[S+2]([O-])([O-])CC")
        ]

        self._test_no_product(reactants, self.reactant_names)
