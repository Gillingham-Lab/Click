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
            with self.subTest(reactants=(amine, acid, product)):
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
            with self.subTest(reactants=(amine, acid)):
                amine = fromSmiles(amine)
                acid = fromSmiles(acid)

                with self.assertRaises(Click.Exceptions.NoProductError):
                    test_product = Click.AmideCoupling(amine=amine, acid=acid).getProduct()

    # This reactions should give multiple possible products
    def test_two_or_more_products(self):
        reactants = [
            ("NCCNC", "OC(C)=O", ["CNCCNC(C)=O", "NCCN(C(C)=O)C"]),
            ("CNCCNC", "OC(C)=O", ["CNCCN(C)C(C)=O", "CNCCN(C)C(C)=O"]),
        ]

        for amine, acid, products in reactants:
            with self.subTest(reactants=(amine, acid, products)):
                amine = fromSmiles(amine)
                acid = fromSmiles(acid)
                products = [toSmiles(fromSmiles(x)) for x in products]

                # getProduct only expects 1 product; must give an exception otherwise
                with self.assertRaises(Click.Exceptions.AmbiguousProductError):
                    test_products = Click.AmideCoupling(amine=amine, acid=acid).getProduct()

                # Use getProducts to get a list of products.
                test_products = Click.AmideCoupling(amine=amine, acid=acid).getProducts()
                found = 0

                for p in test_products:
                    p_smiles = toSmiles(p)

                    self.assertIn(p_smiles, products)
                    found+=1

                self.assertEqual(found, len(test_products), "Not all products were expected.")
                self.assertEqual(len(test_products), len(products), "Mismatch between expected products and actual products.")


if __name__ == '__main__':
    unittest.main()
