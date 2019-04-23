__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.io import _parseMolData, _parseSMILESData, _getSMILESString, _getSDFString
from chembl_beaker.beaker.utils.io import _getSMARTSString, _getXYZ
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords
from chembl_beaker.beaker.utils.chemical_transformation import _sanitize
from chembl_beaker.beaker.utils.io import _molFromSmarts

# ----------------------------------------------------------------------------------------------------------------------


def _canonicalize_smiles(data, computeCoords=False, in_delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True,
                         sanitize=True, out_delimiter=' ', nameHeader='Name', includeHeader=True, isomericSmiles=False,
                         kekuleSmiles=False):
    return _getSMILESString(_parseSMILESData(data, computeCoords=computeCoords, delimiter=in_delimiter,
                                             smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine,
                                             sanitize=sanitize),
                            delimiter=out_delimiter, nameHeader=nameHeader, includeHeader=includeHeader,
                            isomericSmiles=isomericSmiles, kekuleSmiles=kekuleSmiles)

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2smarts(data, sanitize=True, removeHs=True, strictParsing=True, isomericSmiles=False):
    return _getSMARTSString(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
                            isomericSmiles=isomericSmiles)

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2smiles(data, sanitize=True, removeHs=True, strictParsing=True, delimiter=' ', nameHeader='Name',
                 includeHeader=True, isomericSmiles=False, kekuleSmiles=False):
    return _getSMILESString(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
                            delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader,
                            isomericSmiles=isomericSmiles, kekuleSmiles=kekuleSmiles)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2ctab(data, computeCoords=True, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True, sanitize=True):
    return _getSDFString(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
                                          smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine,
                                          sanitize=sanitize))

# ----------------------------------------------------------------------------------------------------------------------


def _smarts2ctab(data, computeCoords=True, delimiter=' ', sanitize=True):
    mols = []
    for line in data.splitlines():
        if not line:
            continue
        for chunk in line.strip().split(delimiter):
            if not chunk:
                continue
            mols.append(_molFromSmarts(chunk))
    if computeCoords:
        _apply(mols, _computeCoords, True)
    if sanitize:
        _apply(mols, _sanitize)
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2xyz(data, computeCoords=True):
    return _getXYZ(_parseMolData(data), computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _inchi2ctab(inchis):
    mols = _apply(inchis.split(),Chem.MolFromInchi)
    _apply(mols, _computeCoords)
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2inchi(data, sanitize=True, removeHs=True, strictParsing=True):
    return '\n'.join(_apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
                            Chem.MolToInchi))

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2inchiKey(data, sanitize=True, removeHs=True, strictParsing=True):
    inchis = _ctab2inchi(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    return _inchi2inchiKey(inchis)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2inchi(data, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True,
                  sanitize=True):
    return '\n'.join(_apply(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
                                             smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine,
                                             sanitize=sanitize),
                            Chem.MolToInchi))

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2inchiKey(data, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True,
                     sanitize=True):
    inchis = _smiles2inchi(data, computeCoords=computeCoords, delimiter=delimiter, smilesColumn=smilesColumn,
                           nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize)
    return _inchi2inchiKey(inchis)

# ----------------------------------------------------------------------------------------------------------------------


def _inchi2inchiKey(inchis):
    return '\n'.join([Chem.InchiToInchiKey(inch) for inch in inchis.split()])

# ----------------------------------------------------------------------------------------------------------------------
