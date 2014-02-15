__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import StringIO
from rdkit.Chem.Draw import SimilarityMaps
from chembl_beaker.beaker.utils.io import _parseSMILESData, _parseMolData

try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
except:
    matplotlib=None

#-----------------------------------------------------------------------------------------------------------------------

def _similarityMap(ms,params):
    if matplotlib is None:
        raise ValueError('matplotlib not useable')

    fp=params.get('fingerprint','morgan')
    if fp=='morgan':
        rad = int(params.get('radius',2))
        fn = lambda x,i:SimilarityMaps.GetMorganFingerprint(x,i,radius=rad)
    elif fp=='tt':
        fn = SimilarityMaps.GetAPFingerprint
    elif fp=='ap':
        fn = SimilarityMaps.GetTTFingerprint

    w = int(params.get('w',100))
    h = int(params.get('h',100))

    fig,maxv = SimilarityMaps.GetSimilarityMapForFingerprint(ms[0],ms[1],fn,size=(w,h))
    sio = StringIO.StringIO()
    pyplot.savefig(sio,format='png',bbox_inches='tight',dpi=100)

    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

def _smiles2SimilarityMap(data,params):
    return _similarityMap(_parseSMILESData(data, True), params)

#-----------------------------------------------------------------------------------------------------------------------

def _sdf2SimilarityMap(data,params):
    return _similarityMap(_parseMolData(data), params)

#-----------------------------------------------------------------------------------------------------------------------
