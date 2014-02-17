__author__ = 'mnowotka'

from bottle import request, response
import json

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.marvin.impl import _clean2D, _stereoInfo, _molExport

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/clean', method=['OPTIONS', 'POST'], name="clean")
def clean2D():
    """
Implements Marvin 4 js clean2D web service. Recomputes 2D coordinates of given compound. The compound is in
Marvin's *.mrv format.
    """

    params = json.loads(request.body.read())
    structure = params['structure']
    response.content_type = 'text/plain'
    return _clean2D(structure)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/cipStereoInfo', method=['OPTIONS', 'POST'], name="cipStereoInfo")
def stereoInfo():
    """
Implements Marvin 4 js Stereo Info web service. Marks possible stereocenters, R/S chirality for atoms and cistrans
bonds.
    """

    params = json.loads(request.body.read())
    structure = params['structure']
    response.content_type = 'application/json'
    return json.dumps(_stereoInfo(structure))

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/molExport' , method=['OPTIONS', 'POST'], name="molExport")
def molExport():
    """
Implements Marvin 4 js MolConvert web service. Converts compound from one format to another, including Marvin's
*.mrv.
    """

    params = json.loads(request.body.read())
    structure = params['structure']
    input_f = params['inputFormat']
    output_f = params['parameters']

    res = _molExport(structure, input_f, output_f)
    response.content_type = 'application/json'
    return json.dumps(res)

#-----------------------------------------------------------------------------------------------------------------------
