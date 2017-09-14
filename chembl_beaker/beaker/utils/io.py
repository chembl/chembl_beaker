__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

import os
import tempfile
import StringIO
from rdkit.Chem import SDMolSupplier
from rdkit.Chem import SmilesMolSupplier
from rdkit.Chem import MolToSmarts
from rdkit.Chem import MolFromSmarts
from rdkit.Chem import SDWriter
from rdkit.Chem import SmilesWriter
from chembl_beaker.beaker.utils.functional import _apply, _call
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords
from chembl_beaker.beaker.utils.chemical_transformation import _getSubstructMatch

# ----------------------------------------------------------------------------------------------------------------------


def _parseFlag(data):
    try:
        return bool(int(data))
    except ValueError:
        if isinstance(data, basestring):
            l = data.lower()
            if l == 't' or l == 'true':
                return True
            if l == 'f' or l == 'false':
                return False
        if not data:
            return False
        return None

# ----------------------------------------------------------------------------------------------------------------------


def _parseMolData(data, sanitize=True, removeHs=True, strictParsing=True):
    fd, fpath = tempfile.mkstemp(text=True)
    os.write(fd, data)
    os.close(fd)
    suppl = SDMolSupplier(fpath, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    res = [x for x in suppl if x]
    os.remove(fpath)
    return res

# ----------------------------------------------------------------------------------------------------------------------


def _parseSMILESData(data, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True,
                     sanitize=True):
    fd, fpath = tempfile.mkstemp(text=True)
    os.write(fd, data)
    os.close(fd)
    suppl = SmilesMolSupplier(fpath, delimiter=delimiter, smilesColumn=smilesColumn, nameColumn=nameColumn,
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
    w = SDWriter(f)
    for m in mols:
        w.write(m)
    w.flush()

# ----------------------------------------------------------------------------------------------------------------------


def _getSDFString(mols):
    sio = StringIO.StringIO()
    _getSDFStream(sio, mols)
    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _getSMILESStream(f, mols, delimiter=' ', nameHeader='Name', includeHeader=True, isomericSmiles=False,
                     kekuleSmiles=False):
    w = SmilesWriter(f, delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader,
                     isomericSmiles=isomericSmiles, kekuleSmiles=kekuleSmiles)
    for mol in mols:
        w.write(mol)
    w.flush()

# ----------------------------------------------------------------------------------------------------------------------


def _getSMILESString(mols, delimiter=' ', nameHeader='Name', includeHeader=True, isomericSmiles=False,
                     kekuleSmiles=False):
    sio = StringIO.StringIO()
    _getSMILESStream(sio, mols, delimiter=delimiter, nameHeader=nameHeader, includeHeader=includeHeader,
                     isomericSmiles=isomericSmiles, kekuleSmiles=kekuleSmiles)
    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _getSMARTSString(mols, isomericSmiles=False):
    sio = StringIO.StringIO()
    for mol in mols:
        sio.write(MolToSmarts(mol, isomericSmiles) + '\n')
    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _molFromSmarts(smarts):
    return MolFromSmarts(str(smarts))

# ----------------------------------------------------------------------------------------------------------------------


def _getXYZ(mols, computeCoords=True):

    sio = StringIO.StringIO()

    for i, mol in enumerate(mols):

        molH = mol
        if computeCoords:
            from chembl_beaker.beaker.core_apps.D3Coords.impl import _2D23D
            molH = _2D23D(mol, None)

        atoms = molH.GetAtoms()
        sio.write(str(len(atoms)) + '\n')
        sio.write('\n')

        i = 0
        for conf in molH.GetConformers():

            for j in range(0, conf.GetNumAtoms()):

                sio.write(atoms[i].GetSymbol() +
                          '\t' + str(conf.GetAtomPosition(j).x) +
                          '\t' + str(conf.GetAtomPosition(j).y) +
                          '\t' + str(conf.GetAtomPosition(j).z)
                          + '\n')
                i += 1

    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _getMatches(mols, smarts, force=False):
    _call(mols, 'UpdatePropertyCache', strict=False)
    return _apply(mols, _getSubstructMatch, _molFromSmarts(smarts), force)

# ----------------------------------------------------------------------------------------------------------------------

