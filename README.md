# Click

Click is a collection of pre-made and tested reaction patterns to use with RDKit molecules.
The reactions can be used to sequentially modify a molecule, or to create combinatorial 
libraries. Click is only doing the reaction as you've specified and does *not* check 
if the reaction would work.

## Requirements

* RDKit (version >= 2019.03)
* Python (version >= 3.6)

## Installation
To install Click, run

`pip install ClickReaction`

## Usage

Many examples can be found in the [tests folder](https://github.com/Gillingham-Lab/Click/tree/master/tests).

### Boc removal

```python
from rdkit import Chem
from ClickReaction import BocRemoval

boc_protected_amine = Chem.MolFromSmiles("CNC(OC(C)(C)C)=O")

reaction = BocRemoval(bocamine=boc_protected_amine)
product = reaction.get_product()

assert "CN" == Chem.MolToSmiles(product)
```

### Click Reaction

```python
from rdkit import Chem
from ClickReaction import CuAAC

alkyne = Chem.MolFromSmiles("c1ccccc1C#C")
azide = Chem.MolFromSmiles("C-[N-]-[N+]#N")

reaction = CuAAC(alkyne=alkyne, azide=azide)
product = reaction.get_product()

assert "Cn1cc(-c2ccccc2)nn1" == Chem.MolToSmiles(product)
```

## Supported reactions

### Simple transformations

* Boc removal
* Fmoc removal
* Alkaline ester hydrolysis

### Bimolecular reactions

* Amide coupling (with or without anilines)
* CuAAC
* Sulfon amide formation from amines and sulfonyl chlorides
* Suzuki-Miyaura cross coupling