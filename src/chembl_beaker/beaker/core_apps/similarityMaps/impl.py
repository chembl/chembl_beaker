__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------


import io
from rdkit.Chem.Draw import SimilarityMaps
from chembl_beaker.beaker.utils.functional import _apply, _call
from chembl_beaker.beaker.utils.io import _parseSMILESData, _parseMolData
import chembl_beaker.beaker.utils.chemical_transformation as ct

try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
except:
    matplotlib=None

# ----------------------------------------------------------------------------------------------------------------------


def _similarityMap(ms, width=500, height=500, radius=2, fingerprint='morgan', format='png'):
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


def _smiles2SimilarityMap(data, width=500, height=500, radius=2, fingerprint='morgan', computeCoords=False,
                          delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True, sanitize=True, format='png'):
    return _similarityMap(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize),
        width=width, height=height, radius=radius, fingerprint=fingerprint, format=format)


# ----------------------------------------------------------------------------------------------------------------------

def _sdf2SimilarityMap(data, width=500, height=500, radius=2, fingerprint='morgan', sanitize=True, removeHs=True,
                       strictParsing=True, format='png'):
    return _similarityMap(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
        width=width, height=height, radius=radius, fingerprint=fingerprint, format=format)

# ----------------------------------------------------------------------------------------------------------------------

