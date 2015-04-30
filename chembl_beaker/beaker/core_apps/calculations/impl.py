__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem.rdmolops import SanitizeFlags as sf
SANITIZE_ALL = sf.SANITIZE_ALL
from chembl_beaker.beaker.utils.functional import _apply, _call
from chembl_beaker.beaker.utils.io import _parseMolData, _getSDFString

#-----------------------------------------------------------------------------------------------------------------------

def _kekulize(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    _apply(mols, Chem.Kekulize)
    return _getSDFString(mols)

#-----------------------------------------------------------------------------------------------------------------------

def _sanitize(data, sanitizeOps=SANITIZE_ALL):
    mols = _parseMolData(data, sanitize=False, removeHs=False, strictParsing=False)
    try:
        _apply(mols, Chem.SanitizeMol, sanitizeOps=sanitizeOps)
    except:
        pass
    return _getSDFString(mols)

#-----------------------------------------------------------------------------------------------------------------------

def _addHs(data, explicitOnly=False, addCoords=False):
    mols = _parseMolData(data, sanitize=True, removeHs=False, strictParsing=True)
    mols = _apply(mols, Chem.AddHs, explicitOnly=explicitOnly, addCoords=addCoords)
    return _getSDFString(mols)

#-----------------------------------------------------------------------------------------------------------------------

def _removeHs(data, implicitOnly=False):
    mols = _parseMolData(data, sanitize=False, removeHs=False, strictParsing=True)
    mols = _apply(mols, Chem.RemoveHs, implicitOnly=implicitOnly)
    return _getSDFString(mols)

#-----------------------------------------------------------------------------------------------------------------------

def _sssr(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return _apply(mols, Chem.GetSSSR)

#-----------------------------------------------------------------------------------------------------------------------

def _symmsssr(data, sanitize=True, removeHs=True, strictParsing=True):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return [[list(xx) for xx in x] for x in _apply(mols, Chem.GetSymmSSSR)]

#-----------------------------------------------------------------------------------------------------------------------