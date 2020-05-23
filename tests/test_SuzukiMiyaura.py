from Click.Reactions import SuzukiMiyaura
from tests import TestHelper


class SulfonAmideFormationTest(TestHelper.ReactionTestCase):
    reactant_names = ["boronate", "halogenide"]

    def setUp(self):
        self.set_reaction(SuzukiMiyaura)

    def test_one_product(self):
        reactants = [
            # Test phenylhalogenides and phenylboronic acid
            (("c1ccccc1B(O)O", "c1ccccc1Cl"), "c1ccccc1c1ccccc1"),
            (("c1ccccc1B(O)O", "c1ccccc1Br"), "c1ccccc1c1ccccc1"),
            (("c1ccccc1B(O)O", "c1ccccc1I"), "c1ccccc1c1ccccc1"),

            # Check different esters for the boronic acid
            # BPin
            (("CC1(C)C(C)(C)OB(C2=CC=CC=C2)O1", "c1ccccc1Cl"), "c1ccccc1c1ccccc1"),
            # Mida (open form)
            (("CN1CC(=O)OB(OC(=O)C1)c2ccccc2", "c1ccccc1Cl"), "c1ccccc1c1ccccc1"),
            # Mida (closed form)
            (("C[N+]12CC(O[B-](c3ccccc3)1OC(C2)=O)=O", "c1ccccc1Cl"), "c1ccccc1c1ccccc1"),

            # Check trifluoroboronates
            (("F[B-](F)(F)C1=CC=CC=C1.[K+]", "c1ccccc1Cl"), "c1ccccc1c1ccccc1"),
        ]

        self._test_one_product(reactants, self.reactant_names)

    # Those tests should not give any product.
    def test_no_product(self):
        reactants = [
            # Arylfluorides should not work.
            ("c1ccccc1B(O)O", "c1ccccc1F"),
            # This should also not work.
            ("c1ccccc1P(O)O", "c1ccccc1Cl"),
        ]

        self._test_no_product(reactants, self.reactant_names)
