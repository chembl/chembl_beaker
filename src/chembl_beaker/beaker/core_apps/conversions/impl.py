__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from beaker.utils.functional import _apply
from beaker.utils.io import _parseMolData, _parseSMILESData, _getSMILESString, _getSDFString
from beaker.utils.io import _getSMARTSString
from beaker.utils.chemical_transformation import _computeCoords
from beaker.utils.chemical_transformation import _sanitize
from beaker.utils.io import _molFromSmarts

# ----------------------------------------------------------------------------------------------------------------------


def _canonicalize_smiles(data, computeCoords=False, in_delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True,
                         sanitize=True, out_delimiter=' ', nameHeader='Name', includeHeader=True):
    return _getSMILESString(_parseSMILESData(data, computeCoords=computeCoords, delimiter=in_delimiter,
                                             smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine,
                                             sanitize=sanitize),
                            delimiter=out_delimiter, nameHeader=nameHeader, includeHeader=includeHeader)

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2smarts(data, loadMol=True, useRDKitChemistry=True):
    return _getSMARTSString(_parseMolData(data, loadMol=True, useRDKitChemistry=True))

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2smiles(data, loadMol=True, useRDKitChemistry=True, delimiter=' ', nameHeader='Name', includeHeader=True):
    return _getSMILESString(_parseMolData(data, loadMol=True, useRDKitChemistry=True),
                            delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2ctab(data, computeCoords=True, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True, sanitize=True):
    return _getSDFString(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
                                          smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine,
                                          sanitize=sanitize))

# ----------------------------------------------------------------------------------------------------------------------


def _smarts2ctab(data, computeCoords=True, delimiter=' ', sanitize=True):
    mols = []
    data = data.decode("utf-8")
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


def _inchi2ctab(inchis):
    mols = _apply(inchis.split(), Chem.MolFromInchi)
    _apply(mols, _computeCoords)
    return _getSDFString(mols)

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2inchi(data, loadMol=False, useRDKitChemistry=False):
    return '\n'.join(_apply(_parseMolData(data, loadMol=False, useRDKitChemistry=False), Chem.MolBlockToInchi))

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2inchiKey(data, loadMol=False, useRDKitChemistry=False):
    inchis = _ctab2inchi(data, loadMol=False, useRDKitChemistry=False)
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
