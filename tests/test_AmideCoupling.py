import Click
from tests import TestHelper


class AmideCouplingTest(TestHelper.ReactionTestCase):
    reactant_names = ["amine", "acid"]

    def setUp(self):
        self.set_reaction(Click.AmideCoupling)

    def test_one_product(self):
        reactants = [
            # Ammonia and ammonium with acetic acid
            (("N", "OC(C)=O"), "O=C(C)N"),
            (("[NH4+]", "OC(C)=O"), "O=C(C)N"),
            (("N(-[H])(-[H])(-[H])", "OC(C)=O"), "O=C(C)N"),

            # Primary amines with simple acids
            (("CN", "OC=O"), "CNC=O"),
            (("CCN", "OC(C)=O"), "O=C(C)NCC"),
            (("CN", "OC(C1=CC=CC=C1)=O"), "O=C(C1=CC=CC=C1)NC"),
            (("NC1=CC=CC(CN)=C1", "OC(C)=O"), "NC1=CC=CC(CNC(C)=O)=C1"),

            # Secondary amines with simple acids
            (("CNC", "OC(C)=O"), "O=C(C)N(C)C"),

            # Primary ammonium with simple acids
            (("CC[NH3+]", "OC(C)=O"), "O=C(C)NCC"),

            # Secondary ammonium with simple acids
            (("C[NH2+]C", "OC(C)=O"), "O=C(C)N(C)C"),

            # Simple amine with carboxylate
            (("CCN", "[O-]C(C)=O"), "O=C(C)NCC"),

            # Secondary ammonium with carboxylate
            (("C[NH2+]C", "[O-]C(C(C1=CC=CC=C1)(F)Cl)=O"), "O=C(C(C1=CC=CC=C1)(F)Cl)N(C)C"),

            #
            # Amine with active esters
            #

            # N-Hydroxysuccinimide
            (("CCN", "CCC(ON1C(CCC1=O)=O)=O"), "CCC(NCC)=O"),
            # N-Sulfo-Hydroxysuccinimide
            (("CCN", "CCC(ON1C(CC(S(=O)(O)=O)C1=O)=O)=O"), "CCC(NCC)=O"),
            # p-Nitrophenol
            (("CCN", "CCC(OC1=CC=C([N+]([O-])=O)C=C1)=O"), "CCC(NCC)=O"),
            # Pentafluorophenyl
            (("CCN", "CCC(OC1=C(F)C(F)=C(F)C(F)=C1F)=O"), "CCC(NCC)=O"),

        ]

        self._test_one_product(reactants, self.reactant_names)

    # Those tests should not give any product.
    def test_no_product(self):
        reactants = [
            # No primary amide
            ("NC(C)=O", "OC(C)=O"),
            # No secondary amide
            ("CC(NC)=O", "OC(C)=O"),
            # No N-methyl urea
            ("NC(NC)=O", "OC(C)=O"),
            # No carbamates
            ("O=C(OC)NC", "OC(C)=O"),
            # No imides
            ("CC(NC(C)=O)=O", "OC(C)=O"),
            # No aniline
            ("NC1=CC=CC=C1", "OC(C1=CC=CC=C1)=O"),
            # No tertiary amines
            ("CN(C)C", "OC(C)=O"),
            # No guanosine
            ("NC(=N)N", "OC(C)=O")
        ]

        self._test_no_product(reactants, self.reactant_names)

    # This reactions should give multiple possible products
    def test_if_get_products_returns_all_possible_products(self):
        reactants = [
            (("NCCNC", "OC(C)=O"), ("CNCCNC(C)=O", "NCCN(C(C)=O)C")),
            (("CNCCNC", "OC(C)=O"), ("CNCCN(C)C(C)=O", "CNCCN(C)C(C)=O")),
            (("N", "OC(CN(CC(O)=O)CCN(CC(O)=O)CC(O)=O)=O"), ["OC(CN(CC(O)=O)CCN(CC(N)=O)CC(O)=O)=O"]*4),
        ]

        self._test_all_possible_products(reactants, self.reactant_names)

    def test_if_get_products_returns_only_1_products_if_symmetrical_as_one_is_true(self):
        reactants = [
            (("NCCNC", "OC(C)=O"), ["CNCCNC(C)=O", "NCCN(C(C)=O)C"]),
            (("CNCCNC", "OC(C)=O"), ["CNCCN(C)C(C)=O"]),
            (("N", "OC(CN(CC(O)=O)CCN(CC(O)=O)CC(O)=O)=O"), ["OC(CN(CC(O)=O)CCN(CC(N)=O)CC(O)=O)=O"]),
        ]

        self._test_all_possible_products(reactants, self.reactant_names, symmetrical_as_one=True)

    def test_if_get_product_returns_the_symmetric_product_if_symmetrical_as_one_is_true(self):
        reactants = [
            (("CNCCNC", "OC(C)=O"), "CNCCN(C)C(C)=O"),
            (("N", "OC(CN(CC(O)=O)CCN(CC(O)=O)CC(O)=O)=O"), "OC(CN(CC(O)=O)CCN(CC(N)=O)CC(O)=O)=O"),
        ]

        self._test_one_product(reactants, self.reactant_names, symmetrical_as_one=True)
