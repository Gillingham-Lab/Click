from Click.Reactions import AlkalineEsterHydrolysis
from tests import TestHelper


class AlkalineEsterHydrolysisTest(TestHelper.ReactionTestCase):
    reactant_names = ["ester"]

    def setUp(self):
        self.set_reaction(AlkalineEsterHydrolysis)

    def test_one_product(self):
        reactants = [
            # Esters from formiates
            (("C(=O)OC", ), "C(=O)O"),
            (("C(=O)OCC", ), "C(=O)O"),

            # acetates
            (("CC(=O)OC", ), "CC(=O)O"),
            (("CC(=O)OCC", ), "CC(=O)O"),
        ]

        self._test_one_product(reactants, self.reactant_names)

    # Those tests should not give any product.
    def test_no_product(self):
        reactants = [
            # tBu should not be removed
            ("C(=O)OC(C)(C)C", ),

            # Bn should not be removed
            ("C(=O)OCc1ccccc1", ),
        ]

        self._test_no_product(reactants, self.reactant_names)
