__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

import chembl_beaker.beaker.utils.chemical_transformation as ct
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _getSDFString

# ----------------------------------------------------------------------------------------------------------------------


def _kekulize(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    _apply(mols, ct._kekulize)
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------


def _sanitize(data, sanitizeOps=ct.SANITIZE_ALL):
    mols = _parseMolData(data, sanitize=False, removeHs=False, strictParsing=False)
    try:
        _apply(mols, ct._sanitize, sanitizeOps=sanitizeOps)
    except:
        pass
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------


def _addHs(data, explicitOnly=False, addCoords=False):
    mols = _parseMolData(data, sanitize=True, removeHs=False, strictParsing=True)
    mols = _apply(mols, ct._addHs, explicitOnly=explicitOnly, addCoords=addCoords)
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------


def _removeHs(data, implicitOnly=False):
    mols = _parseMolData(data, sanitize=False, removeHs=False, strictParsing=True)
    mols = _apply(mols, ct._removeHs, implicitOnly=implicitOnly)
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------


def _sssr(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return _apply(mols, ct._sssr)

# ----------------------------------------------------------------------------------------------------------------------


def _symmsssr(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return [[list(xx) for xx in x] for x in _apply(mols, ct._symmsssr)]

# ----------------------------------------------------------------------------------------------------------------------


def _align(data, template):
    mol1 = _parseMolData(data)
    mol2 = _parseMolData(template)
    pattern = None
    mols = None
    if len(mol1) == 1:
        pattern = mol1[0]
        mols = mol2
    elif len(mol2) == 1:
        pattern = mol2[0]
        mols = mol1
    if not pattern or len(mols) < 2:
        return
    _apply(mols, ct._align, pattern)
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------
