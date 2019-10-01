__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

import os
import tempfile
import io
from rdkit import Chem
from chembl_beaker.beaker.utils.functional import _apply, _call
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords
from chembl_beaker.beaker.utils.chemical_transformation import _getSubstructMatch

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
    """
    Yields molblocks from a sdf file.
    
    :param filename: 
    :return: 
    """
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


def _reapply_molblock_wedging(m):
    for b in m.GetBonds():
        # only do the wedgeing if the bond doesn't already have something there:
        if b.GetBondDir() == Chem.BondDir.NONE and b.HasProp(
                "_MolFileBondStereo"):
            val = b.GetProp("_MolFileBondStereo")
            if val == '1':
                b.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif val == '6':
                b.SetBondDir(Chem.BondDir.BEGINDASH)

# ----------------------------------------------------------------------------------------------------------------------


def _parseMolData(data, sanitize=True, removeHs=True, strictParsing=True, rdkload=True):
    fd, fpath = tempfile.mkstemp(text=True)
    os.write(fd, data)
    os.close(fd)
    suppl = _parse_sdf(fpath)
    res = []
    for molblock in suppl:
        if rdkload:
            mol = Chem.MolFromMolBlock(molblock, sanitize=sanitize, removeHs=sanitize, strictParsing=sanitize)
            if mol:
                _reapply_molblock_wedging(mol)
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
#    if not mols:
#        mols = [MolFromSmiles(data, sanitize=sanitize)]
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


def _getSMILESStream(f, mols, delimiter=' ', nameHeader='Name', includeHeader=True, isomericSmiles=False,
                     kekuleSmiles=False):
    w = Chem.SmilesWriter(f, delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader,
                          isomericSmiles=isomericSmiles, kekuleSmiles=kekuleSmiles)
    for mol in mols:
        w.write(mol)
    w.flush()

# ----------------------------------------------------------------------------------------------------------------------


def _getSMILESString(mols, delimiter=' ', nameHeader='Name', includeHeader=True, isomericSmiles=False,
                     kekuleSmiles=False):
    sio = io.StringIO()
    _getSMILESStream(sio, mols, delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader,
                     isomericSmiles=isomericSmiles, kekuleSmiles=kekuleSmiles)
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
    return _apply(mols, _getSubstructMatch, _molFromSmarts(smarts), force)

# ----------------------------------------------------------------------------------------------------------------------

