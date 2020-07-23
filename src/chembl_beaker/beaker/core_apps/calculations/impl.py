import beaker.utils.chemical_transformation as ct
from beaker.utils.functional import _apply
from beaker.utils.io import _parseMolData, _getSDFString

# ----------------------------------------------------------------------------------------------------------------------

def _addHs(data, explicitOnly=False, addCoords=False):
    mols = _parseMolData(data, loadMol=True, useRDKitChemistry=True)
    mols = _apply(mols, ct._addHs, explicitOnly=explicitOnly, addCoords=addCoords)
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------


def _removeHs(data, implicitOnly=False):
    mols = _parseMolData(data, loadMol=True, useRDKitChemistry=False)
    mols = _apply(mols, ct._removeHs, implicitOnly=implicitOnly)
    return _getSDFString(mols)
