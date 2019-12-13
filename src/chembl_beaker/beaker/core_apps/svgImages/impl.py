__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

import cairo
from itertools import cycle, islice
from rdkit import Chem
from beaker.utils.functional import _apply, _call
from beaker.utils.io import _parseMolData, _parseSMILESData
from beaker.utils.chemical_transformation import _computeCoords, _atomMapNumber, _kekulize
from beaker.utils.io import _getMatches
from rdkit.Chem.Draw import rdMolDraw2D
import beaker.utils.chemical_transformation as ct
from rdkit.Chem.Draw import SimilarityMaps
import io

try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
except:
    matplotlib=None



# ----------------------------------------------------------------------------------------------------------------------


def _mols2svg(mols, size, kekulize=False, atomMapNumber=False, computeCoords=False, highlightAtomLists=None):

    if not mols:
        return ''
    else:
        # process only one molecule
        mols = mols[0:1]

    _call(mols, 'UpdatePropertyCache', strict=False)
    _apply(mols, ct._sssr)

    if computeCoords:
        _apply(mols, _computeCoords, True)
    if atomMapNumber:
        _apply(mols, _atomMapNumber)

    if kekulize:
        _apply(mols, _kekulize)

    if highlightAtomLists:
        highlightAtomLists = highlightAtomLists[0]

    drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    opts = drawer.drawOptions()
    opts.clearBackground=False
    drawer.DrawMolecule(mols[0], highlightAtoms=highlightAtomLists)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2svg(data, size, loadMol=True, useRDKitChemistry=False, kekulize=False,
              atomMapNumber=False, computeCoords=False):
    return _mols2svg(_parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry),
        size, kekulize=kekulize, atomMapNumber=atomMapNumber,
        computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2svg(data, size, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
                titleLine=True, sanitize=True, kekulize=True, atomMapNumber=False):
    return _mols2svg(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize), size,
        kekulize=kekulize, atomMapNumber=atomMapNumber,
        computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _inchi2svg(inchis, size, kekulize=True, atomMapNumber=False,
               computeCoords=False):
    mols = _apply(inchis.split(), Chem.MolFromInchi)
    _apply(mols, _computeCoords)
    return _mols2svg(mols, size, kekulize=kekulize,
                        atomMapNumber=atomMapNumber, computeCoords=computeCoords)

# ----------------------------------------------------------------------------------------------------------------------


def _highlightSmilesFragmentSVG(data, smarts, size, computeCoords=False, delimiter=' ', smilesColumn=0,
                                nameColumn=1, titleLine=True, sanitize=True, kekulize=True,
                                atomMapNumber=False, force=False):
    mols = _parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
                            smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize)
    matches = _getMatches(mols, smarts, force)
    return _mols2svg(mols, size, kekulize=kekulize,
                     atomMapNumber=atomMapNumber, computeCoords=computeCoords, highlightAtomLists=matches)

# ----------------------------------------------------------------------------------------------------------------------


def _highlightCtabFragmentSVG(data, smarts, size, loadMol=True, useRDKitChemistry=False,
                              kekulize=False, atomMapNumber=False, computeCoords=False,
                              force=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    matches = _getMatches(mols, smarts, force)
    return _mols2svg(mols, size, kekulize=kekulize,
                     atomMapNumber=atomMapNumber, computeCoords=computeCoords, highlightAtomLists=matches)

# ----------------------------------------------------------------------------------------------------------------------


def _similarityMapSVG(ms, width=500, height=500, radius=2, fingerprint='morgan', format='svg'):
    if matplotlib is None:
        raise ValueError('matplotlib not useable')

    _call(ms, 'UpdatePropertyCache', strict=False)
    _apply(ms, ct._sssr)

    fn = None
    if fingerprint == 'morgan':
        fn = lambda x, i: SimilarityMaps.GetMorganFingerprint(x, i, radius=radius)
    elif fingerprint == 'tt':
        fn = SimilarityMaps.GetAPFingerprint
    elif fingerprint == 'ap':
        fn = SimilarityMaps.GetTTFingerprint

    SimilarityMaps.GetSimilarityMapForFingerprint(ms[0], ms[1], fn, size=(width, height))
    sio = io.StringIO()
    pyplot.savefig(sio, format=format, bbox_inches='tight', dpi=100)

    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2SimilarityMapSVG(data, width=500, height=500, radius=2, fingerprint='morgan', computeCoords=False,
                          delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True, sanitize=True, format='svg'):
    return _similarityMapSVG(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize),
        width=width, height=height, radius=radius, fingerprint=fingerprint, format=format)


# ----------------------------------------------------------------------------------------------------------------------

def _sdf2SimilarityMapSVG(data, width=500, height=500, radius=2, fingerprint='morgan', loadMol=True, useRDKitChemistry=True, format='svg'):
    return _similarityMapSVG(_parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry),
        width=width, height=height, radius=radius, fingerprint=fingerprint, format=format)

# ----------------------------------------------------------------------------------------------------------------------
