__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

import os
import tempfile
import io
from rdkit import Chem
from beaker.utils.functional import _apply, _call
from beaker.utils.chemical_transformation import _computeCoords
from beaker.utils.chemical_transformation import _getSubstructMatch
import beaker.utils.chemical_transformation as ct
from chembl_structure_pipeline.standardizer import parse_molblock
from bottle import HTTPError

# ----------------------------------------------------------------------------------------------------------------------


def _parseFlag(data):
    try:
        return bool(int(data))
    except ValueError:
        if isinstance(data, str):
            l = data.lower()
            if l == 't' or l == 'true':
                return True
            if l == 'f' or l == 'false':
                return False
        if not data:
            return False
        return None

# ----------------------------------------------------------------------------------------------------------------------


def _parse_sdf(filename):
    with open(filename, "r") as sdf_file:
        molblock = []
        for line in sdf_file:
            line = line.rstrip("\r\n")
            if not line == "$$$$":
                molblock.append(line)
            else:
                if molblock:
                    yield "\n".join(molblock)
                    molblock = []
        # last molecule without "$$$$" line
        if molblock:
            yield "\n".join(molblock)

# ----------------------------------------------------------------------------------------------------------------------


def _parseMolData(data, loadMol=False, useRDKitChemistry=False):
    fd, fpath = tempfile.mkstemp(text=True)
    os.write(fd, data)
    os.close(fd)
    suppl = _parse_sdf(fpath)
    res = []
    for molblock in suppl:
        if loadMol:
            mol = parse_molblock(molblock, useRDKitChemistry=useRDKitChemistry)
            if not mol:
                raise HTTPError(422, "Unprocessable Entity")
            res.append(mol)
        else:
            res.append(molblock)
    os.remove(fpath)
    return res

# ----------------------------------------------------------------------------------------------------------------------


def _parseSMILESData(data, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True,
                     sanitize=True):
    fd, fpath = tempfile.mkstemp(text=True)
    os.write(fd, data)
    os.close(fd)
    suppl = Chem.SmilesMolSupplier(fpath, delimiter=delimiter, smilesColumn=smilesColumn, nameColumn=nameColumn,
                                   titleLine=titleLine, sanitize=sanitize)
    mols = [x for x in suppl if x]
    if computeCoords:
        _apply(mols, _computeCoords, True)
    _call(mols, 'SetProp', '_Name', '')
    os.remove(fpath)
    return mols

# ----------------------------------------------------------------------------------------------------------------------


def _getSDFStream(f, mols):
    w = Chem.SDWriter(f)
    for m in mols:
        w.write(m)
    w.flush()

# ----------------------------------------------------------------------------------------------------------------------


def _getSDFString(mols):
    sio = io.StringIO()
    _getSDFStream(sio, mols)
    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _getSMILESStream(f, mols, delimiter=' ', nameHeader='Name', includeHeader=True):
    w = Chem.SmilesWriter(f, delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader)
    for mol in mols:
        w.write(mol)
    w.flush()

# ----------------------------------------------------------------------------------------------------------------------


def _getSMILESString(mols, delimiter=' ', nameHeader='Name', includeHeader=True):
    sio = io.StringIO()
    _getSMILESStream(sio, mols, delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader)
    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _getSMARTSString(mols, isomericSmiles=False):
    sio = io.StringIO()
    for mol in mols:
        sio.write(Chem.MolToSmarts(mol, isomericSmiles) + '\n')
    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _molFromSmarts(smarts):
    return Chem.MolFromSmarts(smarts)

# ----------------------------------------------------------------------------------------------------------------------


def _getMatches(mols, smarts, force=False):
    _call(mols, 'UpdatePropertyCache', strict=False)
    _apply(mols, ct._sssr)
    return _apply(mols, _getSubstructMatch, _molFromSmarts(smarts), force)

# ----------------------------------------------------------------------------------------------------------------------

