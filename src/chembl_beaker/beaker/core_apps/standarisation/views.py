__author__ = 'efelix'

#-----------------------------------------------------------------------------------------------------------------------

from beaker import app
from bottle import request, response
from beaker.utils.io import _parseFlag
from beaker.core_apps.standarisation.impl import _check, _standardize, _get_parent
import base64

#-----------------------------------------------------------------------------------------------------------------------


def getParentView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', False))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    return _get_parent(data, **kwargs)


#-----------------------------------------------------------------------------------------------------------------------

@app.route('/getParent', method=['OPTIONS', 'POST'], name='getParent')
def get_parent():
    """
Remove salt/solvates.
Examples and documentation: []()
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @unsalt.mol ${BEAKER_ROOT_URL}getParent
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}getParent
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return getParentView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------


def standardizeView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', False))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    return _standardize(data, **kwargs)


#-----------------------------------------------------------------------------------------------------------------------

@app.route('/standardize', method=['OPTIONS', 'POST'], name='standardize')
def standardize():
    """
Get standardized molecule.
Examples and documentation: []()
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @standardize.mol ${BEAKER_ROOT_URL}standardize
    curl -X POST -F "file=@standardize.mol" ${BEAKER_ROOT_URL}standardize
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return standardizeView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------


def checkView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', False))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    return _check(data, **kwargs)


#-----------------------------------------------------------------------------------------------------------------------

@app.route('/check', method=['OPTIONS', 'POST'], name='check')
def check():
    """
Check molecule for issues.
Examples and documentation: []()
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @standardize.mol ${BEAKER_ROOT_URL}check
    curl -X POST -F "file=@standardize.mol" ${BEAKER_ROOT_URL}check

    """
    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return checkView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------