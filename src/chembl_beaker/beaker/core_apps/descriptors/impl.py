__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem import Descriptors
try:
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRotatableBonds
except:
    pass
from chembl_beaker.beaker.utils.functional import _call, _apply
from chembl_beaker.beaker.utils.io import _parseMolData

# ----------------------------------------------------------------------------------------------------------------------


def _desc(mol, name):
    if name and isinstance(name, basestring) and hasattr(Descriptors, name):
        return getattr(Descriptors, name)(mol)

# ----------------------------------------------------------------------------------------------------------------------


def _desc_list(mol, names):
    descriptors = dict()
    for name, fn in Descriptors.descList:
        if not names or name in names:
            descriptors[name] = fn(mol)
    if 'MolecularFormula' not in descriptors:
        descriptors['MolecularFormula'] = CalcMolFormula(mol)
        descriptors['NumRotatableBonds'] = CalcNumRotatableBonds(mol)
    return descriptors

# ----------------------------------------------------------------------------------------------------------------------


def _getNumAtoms(data, sanitize=True, removeHs=True, strictParsing=True):
    return _call(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                               strictParsing=strictParsing), "GetNumAtoms")

# ----------------------------------------------------------------------------------------------------------------------


def _getNumBonds(data, sanitize=True, removeHs=True, strictParsing=True):
    return _call(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                               strictParsing=strictParsing), "GetNumBonds")

# ----------------------------------------------------------------------------------------------------------------------


def _getLogP(data, sanitize=True, removeHs=True, strictParsing=True):
    return _apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                                strictParsing=strictParsing), _desc, 'MolLogP')

# -----------------------------------------------------------------------------------------------------------------------


def _getTPSA(data, sanitize=True, removeHs=True, strictParsing=True):
    return _apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                                strictParsing=strictParsing), _desc, 'TPSA')

# ----------------------------------------------------------------------------------------------------------------------


def _getMolWt(data, sanitize=True, removeHs=True, strictParsing=True):
    return _apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                                strictParsing=strictParsing), _desc, 'MolWt')

# ----------------------------------------------------------------------------------------------------------------------


def _getDescriptors(data, ds, sanitize=True, removeHs=True, strictParsing=True):
    if ds:
        ds = ds.split(',')
    return _apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                                strictParsing=strictParsing), _desc_list, ds)

# ----------------------------------------------------------------------------------------------------------------------
