__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _getSDFString

#-----------------------------------------------------------------------------------------------------------------------

def _kekulize(data):
    mols = _parseMolData(data)
    _apply(mols, Chem.Kekulize)
    return _getSDFString(mols)

#-----------------------------------------------------------------------------------------------------------------------
