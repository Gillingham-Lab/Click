import unittest
from rdkit.Chem import MolFromSmiles as fromSmiles
from rdkit.Chem import MolToSmiles as toSmiles
from typing import Sequence, Tuple, Union, Type

import Click
import Click.Exceptions
import Click.BaseReaction

Reactants = Sequence[str]
ReactantNames = Sequence[str]
NoProductTestCase = Reactants
OneProductTestCase = Tuple[Reactants, str]
MultipleProductTestCase = Tuple[Reactants, Sequence[str]]
AnyProductTestCase = Union[NoProductTestCase, OneProductTestCase, MultipleProductTestCase]


class ReactionTestCase(unittest.TestCase):
    def set_reaction(self, reaction: Type):
        self.reaction = reaction

    def prepare_testcases(self,
                          case: AnyProductTestCase,
                          reactant_names: ReactantNames,
                          with_product: bool = True,
                          ):
        product_expected = None

        if with_product is True:
            # Product is always the last one - we use here to canonicalize given smiles.

            if isinstance(case[1], str):
                product_expected = toSmiles(fromSmiles(case[1]))
            else:
                product_expected = [toSmiles(fromSmiles(x)) for x in case[1]]

            reactants = case[0]
        else:
            reactants = case

        # Convert reactants to objects from smiles
        reactants = [fromSmiles(x) for x in reactants]

        # create a dictonary
        reactants = {reactant_names[x]: reactants[x] for x in range(len(reactants))}

        return product_expected, reactants

    def _test_one_product(self,
                          tests: Sequence[OneProductTestCase],
                          reactant_names: ReactantNames,
                          symmetrical_as_one: bool = False,
                          ):
        """ Tests for giving exactly one expected product. """
        for case in tests:
            with self.subTest(case=case):
                product_expected, reactants = self.prepare_testcases(case, reactant_names)

                # Run the reaction
                product = self.reaction(**reactants).getProduct(symmetrical_as_one=symmetrical_as_one)

                # Test if there was any product
                self.assertIsNotNone(product)

                # Test if the product is what we expected
                product_smiles = toSmiles(product)
                self.assertEqual(product_expected, product_smiles)

    def _test_no_product(self,
                         tests: Sequence[NoProductTestCase],
                         reactant_names: ReactantNames,
                         symmetrical_as_one: bool = False,
                         ):
        """ Tests for failing to give a product. """
        for case in tests:
            with self.subTest(case=case):
                _, reactants = self.prepare_testcases(case, reactant_names, with_product=False)

                with self.assertRaises(Click.Exceptions.NoProductError):
                    test_product = self.reaction(**reactants).getProduct(symmetrical_as_one=symmetrical_as_one)

    def _test_all_possible_products(self,
                                    tests: Sequence[MultipleProductTestCase],
                                    reactant_names: ReactantNames,
                                    symmetrical_as_one: bool = False,
                                    ):
        """ Tests if all possible expected products are covered from one reaction. """
        for case in tests:
            with self.subTest(case=case):
                products_expected, reactants = self.prepare_testcases(case, reactant_names)

                # getProduct only expects 1 product - this must give an Exception
                with self.assertRaises(Click.Exceptions.AmbiguousProductError):
                    product = self.reaction(**reactants).getProduct()

                # getProducts should be used instead.
                products = self.reaction(**reactants).getProducts(symmetrical_as_one=symmetrical_as_one)

                products_found = 0
                for p in products:
                    # Convert product to smiles
                    p_smiles = toSmiles(p)

                    self.assertIn(p_smiles, products_expected)
                    products_found += 1

                self.assertEqual(products_found, len(products), "Not all products wer expected")
                self.assertEqual(len(products), len(products_expected), "Mismatch between expected products and actual products.")
