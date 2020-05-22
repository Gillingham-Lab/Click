import Click
from tests import TestHelper


class CuAACTest(TestHelper.ReactionTestCase):
    reactant_names = ["alkyne", "azide"]

    def setUp(self):
        self.set_reaction(Click.CuAAC)

    def test_one_product(self):
        reactants = [
            # S
            (("CC#C", "C-N=[N+]=[N-]"), "Cn1cc(C)nn1"),
            (("CC#C", "C-[N-]-[N+]#N"), "Cn1cc(C)nn1"),
            (("CC#CI", "C-[N-]-[N+]#N"), "Cn1c(I)c(C)nn1"),
            (("c1ccccc1C#C", "C-[N-]-[N+]#N"), "Cn1cc(-c2ccccc2)nn1"),
        ]

        self._test_one_product(reactants, self.reactant_names)

    def test_no_product(self):
        reactants = [
            # Non-terminal alkynes
            ("CC#CC", "C-N=[N+]=[N-]"),
            ("CC#Cc1ccccc1", "C-N=[N+]=[N-]"),
            ("c1ccccc1C#Cc1ccccc1", "C-N=[N+]=[N-]"),
            ("C1CCCC#CCC1", "C-N=[N+]=[N-]"),
        ]

        self._test_no_product(reactants, self.reactant_names)

    def test_if_get_products_returns_all_possible_products(self):
        reactants = [
            (("C#CC#C", "C-N=[N+]=[N-]"), ("C#Cc1cn(C)nn1", "C#Cc1cn(C)nn1")),
            (("C#CCC#CI", "C-N=[N+]=[N-]"), ("IC#CCc1cn(C)nn1", "C#CCC(N=NN1C)=C1I")),
        ]

        self._test_all_possible_products(reactants, self.reactant_names)

    def test_if_get_products_returns_only_1_products_if_symmetrical_as_one_is_true(self):
        reactants = [
            (("C#C", "C-N=[N+]=[N-]"), ["Cn1ccnn1"]),
            (("C#C", "C-[N-]-[N+]#N"), ["Cn1ccnn1"]),
        ]

        self._test_all_possible_products(reactants, self.reactant_names, symmetrical_as_one=True)

    def test_if_get_product_returns_the_symmetric_product_if_symmetrical_as_one_is_true(self):
        reactants = [
            (("C#C", "C-N=[N+]=[N-]"), "Cn1ccnn1"),
            (("C#C", "C-[N-]-[N+]#N"), "Cn1ccnn1"),
        ]

        self._test_one_product(reactants, self.reactant_names, symmetrical_as_one=True)


if __name__ == '__main__':
    unittest.main()
