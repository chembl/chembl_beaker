__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

import os
import tempfile
import io
from rdkit import Chem
from beaker.utils.functional import _apply, _call
from beaker.utils.chemical_transformation import _computeCoords
from beaker.utils.chemical_transformation import _getSubstructMatch

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


def _parse_molblock(molblock, useRDKitChemistry=False):
    m = Chem.MolFromMolBlock(molblock, sanitize=useRDKitChemistry, 
                             removeHs=useRDKitChemistry)
    if not useRDKitChemistry:
        # the RDKit has, by default, removed bond wedging information from the molecule
        # put that back in:
        _reapply_molblock_wedging(m)
        # Set the stereochemistry of double bonds
        # This block can be removed if github #X ends up being accepted and fixed
        anybonds = []
        for bond in m.GetBonds():
            if bond.GetStereo() == Chem.BondStereo.STEREOANY:
                anybonds.append(bond.GetIdx())
        Chem.SetBondStereoFromDirections(m)
        for bidx in anybonds:
            m.GetBondWithIdx(bidx).SetStereo(Chem.BondStereo.STEREOANY)    
    return m

# ----------------------------------------------------------------------------------------------------------------------


def _parseMolData(data, loadMol=False, useRDKitChemistry=False):
    fd, fpath = tempfile.mkstemp(text=True)
    os.write(fd, data)
    os.close(fd)
    suppl = _parse_sdf(fpath)
    res = []
    for molblock in suppl:
        if loadMol:
            mol = _parse_molblock(molblock, useRDKitChemistry=useRDKitChemistry)
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
    return _apply(mols, _getSubstructMatch, _molFromSmarts(smarts), force)

# ----------------------------------------------------------------------------------------------------------------------

