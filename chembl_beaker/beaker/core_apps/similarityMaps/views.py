__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker.utils.io import _parseFlag
from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.similarityMaps.impl import _smiles2SimilarityMap, _sdf2SimilarityMap
from bottle import request, response
import base64

# ----------------------------------------------------------------------------------------------------------------------


def smiles2SimilarityMapView(data, params):

    kwargs = dict()
    kwargs['computeCoords'] = _parseFlag(params.get('computeCoords', False))
    kwargs['delimiter'] = params.get('delimiter', ' ')
    kwargs['smilesColumn'] = int(params.get('smilesColumn', 0))
    kwargs['nameColumn'] = int(params.get('nameColumn', 1))
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['width'] = int(params.get('width', 100))
    kwargs['height'] = int(params.get('height', 100))
    kwargs['radius'] = int(params.get('radius', 2))
    kwargs['fingerprint'] = params.get('fingerprint', 'morgan')
    kwargs['format'] = params.get('format', 'png')

    if params.get('titleLine') is None and not data.startswith('SMILES Name'):
        kwargs['titleLine'] = False
    else:
        kwargs['titleLine'] = _parseFlag(params.get('titleLine', True))

    response.content_type = 'image/png'
    ret = _smiles2SimilarityMap(data, **kwargs)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2SimilarityMap/<smiles>', method=['OPTIONS', 'GET'], name="smiles2SimilarityMap")
def smiles2SimilarityMap(smiles):
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires *exactly two SMILES*
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}smiles2SimilarityMap/$(cat sim.smi | base64 -w 0 | tr "+/" "-_") > sim.png
    curl -X GET "${BEAKER_ROOT_URL}smiles2SimilarityMap/"$(cat sim.smi | base64 -w 0 | tr "+/" "-_")"?width=500&height=500" > sim.png
    curl -X GET "${BEAKER_ROOT_URL}smiles2SimilarityMap/"$(cat sim.smi | base64 -w 0 | tr "+/" "-_")"?width=500&height=500&fingerprint=tt" > sim.png
    curl -X GET "${BEAKER_ROOT_URL}smiles2SimilarityMap/"$(cat sim.smi | base64 -w 0 | tr "+/" "-_")"?width=500&height=500&fingerprint=ap" > sim.png

    """

    data = base64.urlsafe_b64decode(smiles)
    return smiles2SimilarityMapView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/smiles2SimilarityMap', method=['OPTIONS', 'POST'], name="smiles2SimilarityMap")
def smiles2SimilarityMap():
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires *exactly two SMILES*
cURL examples:

    curl -X POST --data-binary @sim.smi ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.png
    curl -X POST -F "file=@sim.smi" -F "width=500" -F "height=500" -F "fingerprint=tt" ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.png
    curl -X POST -F "file=@sim.smi" -F "width=500" -F "height=500" -F "fingerprint=ap" ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.png
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return smiles2SimilarityMapView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def sdf2SimilarityMapView(data, params):

    kwargs = dict()
    kwargs['sanitize'] = _parseFlag(params.get('sanitize', True))
    kwargs['removeHs'] = _parseFlag(params.get('removeHs', True))
    kwargs['strictParsing'] = _parseFlag(params.get('strictParsing', True))
    kwargs['width'] = int(params.get('width', 100))
    kwargs['height'] = int(params.get('height', 100))
    kwargs['radius'] = int(params.get('radius', 2))
    kwargs['fingerprint'] = params.get('fingerprint', 'morgan')
    kwargs['format'] = params.get('format', 'png')

    response.content_type = 'image/png'
    ret = _sdf2SimilarityMap(data, **kwargs)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/sdf2SimilarityMap/<ctab>', method=['OPTIONS', 'GET'], name="sdf2SimilarityMap")
def sdf2SimilarityMap(ctab):
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires SDF containing *exactly two mols*
cURL examples:

    curl -X GET ${BEAKER_ROOT_URL}sdf2SimilarityMap/$(cat sim.sdf | base64 -w 0 | tr "+/" "-_") > sim.png
    curl -X GET "${BEAKER_ROOT_URL}sdf2SimilarityMap/"$(cat sim.sdf | base64 -w 0 | tr "+/" "-_")"?width=500&height=500" > sim.png
    curl -X GET "${BEAKER_ROOT_URL}sdf2SimilarityMap/"$(cat sim.sdf | base64 -w 0 | tr "+/" "-_")"?width=500&height=500&fingerprint=tt" > sim.png
    curl -X GET "${BEAKER_ROOT_URL}sdf2SimilarityMap/"$(cat sim.sdf | base64 -w 0 | tr "+/" "-_")"?width=500&height=500&fingerprint=ap" > sim.png
    """

    data = base64.urlsafe_b64decode(ctab)
    return sdf2SimilarityMapView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/sdf2SimilarityMap', method=['OPTIONS', 'POST'], name="sdf2SimilarityMap")
def sdf2SimilarityMap():
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires SDF containing *exactly two mols*
cURL examples:

    curl -X POST --data-binary @sim.sdf ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.png
    curl -X POST -F "file=@sim.sdf" -F "width=500" -F "height=500" -F "fingerprint=tt" ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.png
    curl -X POST -F "file=@sim.sdf" -F "width=500" -F "height=500" -F "fingerprint=ap" ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.png
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return sdf2SimilarityMapView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------

