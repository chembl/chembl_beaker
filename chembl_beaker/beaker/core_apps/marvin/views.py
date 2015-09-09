__author__ = 'mnowotka'

from bottle import request, response
import json

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.marvin.impl import _clean, _stereoInfo, _molExport, _hydrogenize

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/clean', method=['OPTIONS', 'POST'], name="clean")
def clean():
    """
Implements Marvin 4 js clean2D/3D web service. Recomputes 2D/3D coordinates of given compound. The compound is in
Marvin's *.mrv format. Dim is an optional parameter specifying if 2D or 3D coordinates should be computed
    """

    params = json.loads(request.body.read())
    structure = params['structure']
    dim = params.get('dim', 2)
    response.content_type = 'text/plain'
    return _clean(structure, dim)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/hydrogenize', method=['OPTIONS', 'POST'], name="hydrogenize")
def hydrogenize():
    """
Implements Marvin 4 js Hydrogenize web service. Adds or removes hydrogen atoms. The compound is in
Marvin's *.mrv format. Dim is an optional parameter specifying if 2D or 3D coordinates should be computed
    """

    params = json.loads(request.body.read())
    structure = params['structure']
    input_format = params.get('inputFormat', 'mrv')
    parameters = params.get('parameters')
    method = 'hydrogenize'
    if parameters:
        method = parameters.get('method', 'hydrogenize')
        method = method.lower()
    response.content_type = 'text/plain'
    return _hydrogenize(_molExport(structure, input=input_format, output='mol')["structure"], method == 'hydrogenize')

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
    input_f = params.get('inputFormat')
    output_f = params.get('parameters','mrv')

    res = _molExport(structure, input=input_f, output=output_f)
    response.content_type = 'application/json'
    return json.dumps(res)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/reactionExport' , method=['OPTIONS', 'POST'], name="reactionExport")
def reactionExport():
    """
Not implemented yet
    """
    response.status = 501

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/reactionConverter' , method=['OPTIONS', 'POST'], name="reactionConverter")
def reactionConverter():
    """
Not implemented yet
    """

    response.status = 501

#-----------------------------------------------------------------------------------------------------------------------
