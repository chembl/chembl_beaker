__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from chembl_beaker.beaker.core_apps.fingerprints.impl import _sdf2fps
from bottle import request
import base64
import json

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sdf2fps/<ctab>', method=['OPTIONS', 'GET'], name="sdf2fps")
@app.route('/sdf2fps/<ctab>/<type>', method=['OPTIONS', 'GET'], name="sdf2fps")
def sdf2fps(ctab, type='morgan'):
    """
Computes fingerprints for given compounds and writes them as fps format. Type is a optional type of fingerprints,
can be 'morgan', 'pair' or 'maccs' (default to 'morgan'). CTAB is urlsafe_base64 encoded string containing single
molfile or concatenation of multiple molfiles.
    """
    sdf = base64.urlsafe_b64decode(ctab)
    return _sdf2fps(sdf, type)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sdf2fps', method=['OPTIONS', 'POST'], name="sdf2fps")
def sdf2fps():
    """
Computes fingerprints for given compounds and writes them as fps format. Type is a optional type of fingerprints,
can be 'morgan', 'pair' or 'maccs' (default to 'morgan'). CTAB is either single molfile or SDF file.
    """

    params = json.loads(request.body.read())
    structure = params['structure']
    type = params.get('type', 'morgan')
    radius = params.get('radius', 2)
    n_bits = params.get('n_bits', 2048)
    return _sdf2fps(structure, type, radius, n_bits)

#-----------------------------------------------------------------------------------------------------------------------
