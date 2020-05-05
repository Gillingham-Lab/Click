import os
import sys
from rdkit import Chem

import Click


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def toSmiles(molecule: "rdkit.Chem.rdchem.Mol"):
    return Chem.MolToSmiles(molecule)


def fromSmiles(smiles: str):
    return Chem.MolFromSmiles(smiles)