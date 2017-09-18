__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------


import StringIO
from rdkit.Chem.Draw import SimilarityMaps
from chembl_beaker.beaker.utils.io import _parseSMILESData, _parseMolData

try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
except:
    matplotlib=None

# ----------------------------------------------------------------------------------------------------------------------


def _similarityMap(ms, width=100, height=100, radius=2, fingerprint='morgan', format='png'):
    if matplotlib is None:
        raise ValueError('matplotlib not useable')

    fn = None
    if fingerprint == 'morgan':
        fn = lambda x, i: SimilarityMaps.GetMorganFingerprint(x, i, radius=radius)
    elif fingerprint == 'tt':
        fn = SimilarityMaps.GetAPFingerprint
    elif fingerprint == 'ap':
        fn = SimilarityMaps.GetTTFingerprint

    fig, maxv = SimilarityMaps.GetSimilarityMapForFingerprint(ms[0], ms[1], fn, size=(width, height))
    sio = StringIO.StringIO()
    pyplot.savefig(sio, format=format, bbox_inches='tight', dpi=100)

    return sio.getvalue()

# ----------------------------------------------------------------------------------------------------------------------


def _smiles2SimilarityMap(data, width=100, height=100, radius=2, fingerprint='morgan', computeCoords=False,
                          delimiter=' ', smilesColumn=0, nameColumn=1, titleLine=True, sanitize=True, format='png'):
    return _similarityMap(_parseSMILESData(data, computeCoords=computeCoords, delimiter=delimiter,
        smilesColumn=smilesColumn, nameColumn=nameColumn, titleLine=titleLine, sanitize=sanitize),
        width=width, height=height, radius=radius, fingerprint=fingerprint, format=format)


# ----------------------------------------------------------------------------------------------------------------------

def _sdf2SimilarityMap(data, width=100, height=100, radius=2, fingerprint='morgan', sanitize=True, removeHs=True,
                       strictParsing=True, format='png'):
    return _similarityMap(_parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
        width=width, height=height, radius=radius, fingerprint=fingerprint, format=format)

# ----------------------------------------------------------------------------------------------------------------------

