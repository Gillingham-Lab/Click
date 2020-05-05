import unittest
from rdkit.Chem import MolFromSmiles as fromSmiles
from rdkit.Chem import MolToSmiles as toSmiles

import Click

class AmideCouplingTest(unittest.TestCase):
    def test_one_product(self):
        reactants = [
            ("CN", "OC=O", "CNC=O"),
            ("CCN", "OC(C)=O", "O=C(C)NCC"),
            ("CN", "OC(C1=CC=CC=C1)=O", "O=C(C1=CC=CC=C1)NC"),
        ]

        for amine, acid, product in reactants:
            amine = fromSmiles(amine)
            acid = fromSmiles(acid)
            product = toSmiles(fromSmiles(product))

            test_product = Click.AmideCoupling(amine=amine, acid=acid).getProduct()
            self.assertIsNotNone(test_product)

            test_product_smiles = toSmiles(test_product)
            self.assertEqual(product, test_product_smiles)


if __name__ == '__main__':
    unittest.main()
