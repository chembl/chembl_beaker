__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from itertools import cycle, islice
from chembl_beaker.beaker.utils.functional import _apply
from chembl_beaker.beaker.utils.chemical_transformation import _computeCoords, _atomMapNumber, _kekulize
from chembl_beaker.beaker.utils.io import _parseMolData, _parseSMILESData
from chembl_beaker.beaker.utils.io import _getMatches

NEW_RENDER_ENGINE = False

try:
    from rdkit.Chem.Draw import rdMolDraw2D
    NEW_RENDER_ENGINE = True
except:
    pass
import StringIO
from chembl_beaker.beaker import draw

# ----------------------------------------------------------------------------------------------------------------------


def _mols2imageStream(mols, f, format, size, legend, highlightAtomLists=None, kekulize=True):
    labels = [x for x in islice(cycle(legend), len(mols))] if isinstance(legend, (list, tuple)) else \
             [x for x in islice(cycle([legend]), len(mols))]
    legends = [x.GetProp("_Name") if (x.HasProp("_Name") and x.GetProp("_Name")) else labels[idx]
               for idx, x in enumerate(mols)]

    molsPerRow = min(len(mols), 4)
    nRows = len(mols) // molsPerRow
    if len(mols) % molsPerRow:
        nRows += 1

    if NEW_RENDER_ENGINE:
        if kekulize:
            _apply(mols, _kekulize)

        highlightBondLists = []
        if highlightAtomLists:
            for mol, highlightAtomList in zip(mols, highlightAtomLists):
                highlightBondList = []
                for bnd in mol.GetBonds():
                    if bnd.GetBeginAtomIdx() in highlightAtomList and bnd.GetEndAtomIdx() in highlightAtomList:
                        highlightBondList.append(bnd.GetIdx())
                highlightBondLists.append(highlightBondList)
        if not highlightBondLists:
            highlightBondLists = None

        panelx = size
        panely = size
        canvasx = panelx * molsPerRow
        canvasy = panely * nRows
        drawer = rdMolDraw2D.MolDraw2DCairo(canvasx, canvasy, panelx, panely)
        drawer.DrawMolecules(mols, highlightAtoms=highlightAtomLists, highlightBonds=highlightBondLists, legends=legends)
        drawer.FinishDrawing()
        f.write(drawer.GetDrawingText())

    else:
        image = draw.MolsToGridImage(mols, molsPerRow=molsPerRow, subImgSize=(size, size), legends=legends,
                                     highlightAtomLists=highlightAtomLists)
        image.save(f, format)

# ----------------------------------------------------------------------------------------------------------------------


def _mols2imageString(mols, size, legend, format, atomMapNumber=False, computeCoords=False, highlightAtomLists=None,
                      kekulize=True):
    if not mols:
        return ''
    if computeCoords:
        _apply(mols, _computeCoords, True)
    if atomMapNumber:
        _apply(mols, _atomMapNumber)
    image_data = StringIO.StringIO()
    _mols2imageStream(mols, image_data, format, size, legend, highlightAtomLists, kekulize)
    return image_data.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _ctab2image(data, size, legend, sanitize=True, removeHs=True, strictParsing=True, atomMapNumber=False,
                computeCoords=False, kekulize=True):
    return _mols2imageString(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
                             size, legend, 'PNG', atomMapNumber, computeCoords, kekulize=kekulize)

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2image(data, size, legend, computeCoords=False, delimiter=' ', smilesColumn=0, nameColumn=1,
                  titleLine=True, sanitize=True, atomMapNumber=False, kekulize=True):
    return _mols2imageString(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
                                              smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine,
                                              sanitize=sanitize),
                             size, legend, 'PNG', atomMapNumber, computeCoords, kekulize=kekulize)

# ----------------------------------------------------------------------------------------------------------------------


def _highlightSmilesFragment(data, smarts, size, legend, computeCoords=False, delimiter=' ', smilesColumn=0, 
                             nameColumn=1, titleLine=True, sanitize=True, atomMapNumber=False, kekulize=True,
                             force=False):
    mols = _parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
                            smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize)
    matches = _getMatches(mols, smarts, force)
    return _mols2imageString(mols, size, legend, 'PNG', atomMapNumber, computeCoords, matches, kekulize=kekulize)

# ----------------------------------------------------------------------------------------------------------------------


def _highlightCtabFragment(data, smarts, size, legend, sanitize=True, removeHs=True, strictParsing=True,
                           atomMapNumber=False, computeCoords=False, kekulize=True, force=False):
    mols = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    matches = _getMatches(mols, smarts, force)
    return _mols2imageString(mols, size, legend, 'PNG', atomMapNumber, computeCoords, matches, kekulize=kekulize)

# ----------------------------------------------------------------------------------------------------------------------

