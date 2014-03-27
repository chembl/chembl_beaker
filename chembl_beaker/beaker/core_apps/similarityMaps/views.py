__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.similarityMaps.impl import _smiles2SimilarityMap, _sdf2SimilarityMap
from bottle import request, response
import base64

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2SimilarityMap/<smiles>', method=['OPTIONS', 'GET'], name="smiles2SimilarityMap")
def smiles2SimilarityMap(smiles):
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires *exactly two SMILES*
    """

    data = base64.urlsafe_b64decode(smiles)
    response.content_type = 'image/png'
    ret = _smiles2SimilarityMap(data, request.params)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/smiles2SimilarityMap', method=['OPTIONS', 'POST'], name="smiles2SimilarityMap")
def smiles2SimilarityMap():
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires *exactly two SMILES*
    """

    data = request.body.read()
    response.content_type = 'image/png'
    ret = _smiles2SimilarityMap(data, request.params)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sdf2SimilarityMap/<ctab>', method=['OPTIONS', 'GET'], name="sdf2SimilarityMap")
def sdf2SimilarityMap(ctab):
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires SDF containing *exactly two mols*
    """

    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'image/png'
    ret = _sdf2SimilarityMap(data, request.params)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sdf2SimilarityMap', method=['OPTIONS', 'POST'], name="sdf2SimilarityMap")
def sdf2SimilarityMap():
    """
Generates similarity map, which is a way to visualize the atomic contributions to the similarity between molecules.
This method requires SDF containing *exactly two mols*
    """

    data = request.body.read()
    response.content_type = 'image/png'
    ret = _sdf2SimilarityMap(data, request.params)
    if request.is_ajax:
        ret = base64.b64encode(ret)
    return ret

#-----------------------------------------------------------------------------------------------------------------------
