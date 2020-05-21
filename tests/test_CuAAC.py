import unittest
from rdkit.Chem import MolFromSmiles as fromSmiles
from rdkit.Chem import MolToSmiles as toSmiles

from Click import CuAAC
import Click.Exceptions


class CuAACTest(unittest.TestCase):
    def test_one_product(self):
        reactants = [
            # S
            ("CC#C", "C-N=[N+]=[N-]", "Cn1cc(C)nn1"),
            ("CC#C", "C-[N-]-[N+]#N", "Cn1cc(C)nn1"),
            ("CC#C[I]", "C-[N-]-[N+]#N", "Cn1c([I])c(C)nn1"),
            ("c1ccccc1C#C", "C-[N-]-[N+]#N", "Cn1cc(-c2ccccc2)nn1"),
        ]

        for alkyne, azide, product in reactants:
            with self.subTest(reactants=(alkyne, azide, product)):
                alkyne = fromSmiles(alkyne)
                azide = fromSmiles(azide)
                product = toSmiles(fromSmiles(product))

                test_product = CuAAC(alkyne=alkyne, azide=azide).getProduct()
                self.assertIsNotNone(test_product)

                test_product_smiles = toSmiles(test_product)
                self.assertEqual(product, test_product_smiles)

    # Those tests should not give any product.
    def test_no_product(self):
        reactants = [
            # Non-terminal alkynes
            ("CC#CC", "C-N=[N+]=[N-]"),
            ("CC#Cc1ccccc1", "C-N=[N+]=[N-]"),
            ("c1ccccc1C#Cc1ccccc1", "C-N=[N+]=[N-]"),
        ]

        for alkyne, azide in reactants:
            with self.subTest(reactants=(alkyne, azide)):
                alkyne = fromSmiles(alkyne)
                azide = fromSmiles(azide)

                with self.assertRaises(Click.Exceptions.NoProductError):
                    test_product = CuAAC(alkyne=alkyne, azide=azide).getProduct()

    # This reactions should give multiple possible products
    def test_if_get_products_returns_all_possible_products(self):
        reactants = [
            ("C#CC#C", "C-N=[N+]=[N-]", ["C#Cc1cn(C)nn1", "C#Cc1cn(C)nn1"]),
            ("C#CCC#CI", "C-N=[N+]=[N-]", ["IC#CCc1cn(C)nn1", "C#CCC(N=NN1C)=C1I"])
        ]

        for alkyne, azide, products in reactants:
            with self.subTest(reactants=(alkyne, azide, products)):
                alkyne = fromSmiles(alkyne)
                azide = fromSmiles(azide)
                products = [toSmiles(fromSmiles(x)) for x in products]

                # getProduct only expects 1 product; must give an exception otherwise
                with self.assertRaises(Click.Exceptions.AmbiguousProductError):
                    test_products = CuAAC(alkyne=alkyne, azide=azide).getProduct()

                # Use getProducts to get a list of products.
                test_products = CuAAC(alkyne=alkyne, azide=azide).getProducts()
                found = 0

                for p in test_products:
                    p_smiles = toSmiles(p)

                    self.assertIn(p_smiles, products)
                    found+=1

                self.assertEqual(found, len(test_products), "Not all products were expected.")
                self.assertEqual(len(test_products), len(products), "Mismatch between expected products and actual products.")

    def test_if_get_products_returns_only_1_products_if_symmetrical_as_one_is_true(self):
        reactants = [
            ("C#C", "C-N=[N+]=[N-]", ["Cn1ccnn1"]),
            ("C#C", "C-[N-]-[N+]#N", ["Cn1ccnn1"]),
        ]

        for alkyne, azide, products in reactants:
            with self.subTest(reactants=(alkyne, azide, products)):
                alkyne = fromSmiles(alkyne)
                azide = fromSmiles(azide)
                products = [toSmiles(fromSmiles(x)) for x in products]

                # Use getProducts to get a list of products.
                test_products = CuAAC(alkyne=alkyne, azide=azide).getProducts(symmetrical_as_one=True)
                found = 0

                for p in test_products:
                    p_smiles = toSmiles(p)

                    self.assertIn(p_smiles, products)
                    found+=1

                self.assertEqual(found, len(test_products), "Not all products were expected.")
                self.assertEqual(len(test_products), len(products), "Mismatch between expected products and actual products.")

    def test_if_get_product_returns_the_symmetric_product_if_symmetrical_as_one_is_true(self):
        reactants = [
            ("C#C", "C-N=[N+]=[N-]", "Cn1ccnn1"),
            ("C#C", "C-[N-]-[N+]#N", "Cn1ccnn1"),
        ]

        for alkyne, azide, product in reactants:
            with self.subTest(reactants=(alkyne, azide, product)):
                alkyne = fromSmiles(alkyne)
                azide = fromSmiles(azide)
                product = toSmiles(fromSmiles(product))

                test_product = CuAAC(alkyne=alkyne, azide=azide).getProduct(symmetrical_as_one=True)
                self.assertIsNotNone(test_product)

                test_product_smiles = toSmiles(test_product)
                self.assertEqual(product, test_product_smiles)


if __name__ == '__main__':
    unittest.main()
