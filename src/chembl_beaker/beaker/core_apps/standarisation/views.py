__author__ = 'efelix'

#-----------------------------------------------------------------------------------------------------------------------

from chembl_beaker.beaker import app
from bottle import request, response
from chembl_beaker.beaker.utils.io import _parseFlag
from chembl_beaker.beaker.core_apps.standarisation.impl import _check, _standardise, _get_parent
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

    curl -X POST --data-binary @unsalt.mol ${BEAKER_ROOT_URL}unsalt
    curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}unsalt
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return getParentView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------


def standardiseView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', False))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    return _standardise(data, **kwargs)


#-----------------------------------------------------------------------------------------------------------------------

@app.route('/standardise', method=['OPTIONS', 'POST'], name='standardise')
def standardise():
    """
    Get standardised molecule.
    Examples and documentation: []()
    CTAB is either single molfile or SDF file.
    cURL examples:

    curl -X POST --data-binary @standardise.mol ${BEAKER_ROOT_URL}standardise
    curl -X POST -F "file=@standardise.mol" ${BEAKER_ROOT_URL}standardise
    """

    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return standardiseView(data, request.params)

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
    """
    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return checkView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------