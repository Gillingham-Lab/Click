import unittest
from rdkit.Chem import MolFromSmiles as fromSmiles
from rdkit.Chem import MolToSmiles as toSmiles

import Click
import Click.Exceptions


class AmideCouplingTest(unittest.TestCase):
    def test_one_product(self):
        reactants = [
            # Primary amines
            ("CN", "OC=O", "CNC=O"),
            ("CCN", "OC(C)=O", "O=C(C)NCC"),
            ("N", "OC(C)=O", "O=C(C)N"),
            ("CN", "OC(C1=CC=CC=C1)=O", "O=C(C1=CC=CC=C1)NC"),
            ("NC1=CC=CC(CN)=C1", "OC(C)=O", "NC1=CC=CC(CNC(C)=O)=C1"),

            # Secondary amines
            ("CNC", "OC(C)=O", "O=C(C)N(C)C"),
        ]

        for amine, acid, product in reactants:
            amine = fromSmiles(amine)
            acid = fromSmiles(acid)
            product = toSmiles(fromSmiles(product))

            test_product = Click.AmideCoupling(amine=amine, acid=acid).getProduct()
            self.assertIsNotNone(test_product)

            test_product_smiles = toSmiles(test_product)
            self.assertEqual(product, test_product_smiles)

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
        ]

        for amine, acid in reactants:
            amine = fromSmiles(amine)
            acid = fromSmiles(acid)

            with self.assertRaises(Click.Exceptions.NoProductError):
                test_product = Click.AmideCoupling(amine=amine, acid=acid).getProduct()


if __name__ == '__main__':
    unittest.main()
