
from .Reaction import Reaction

class CuAAC(Reaction):
    """
    Copper-catalyzed alkyne azide cycloaddition to give the 1,4 regioisomer.

    Alkyne + Azide -> 1,4-triazole

    Attributes
    ----------
    alkyne: rdkit.Chem.rdchem.Mol
        A terminal alkyne
    azide: rdkit.Chem.rdchem.Mol
        An azide
    """