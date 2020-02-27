__author__ = 'efelix'

#-----------------------------------------------------------------------------------------------------------------------

from beaker import app
from bottle import request, response
from beaker.utils.io import _parseFlag
from beaker.core_apps.standarisation.impl import _check, _standardize, _get_parent, _exclude

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
Removes salt/solvates using ChEMBL Structure Pipeline functionality.
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
Standardises a molecule using ChEMBL Structure Pipeline functionality.
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
Check molecule for issues using ChEMBL Structure Pipeline functionality.
Examples and documentation: []()
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @check.mol ${BEAKER_ROOT_URL}check
    curl -X POST -F "file=@check.mol" ${BEAKER_ROOT_URL}check

    """
    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return checkView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------


def excludeView(data, params):

    kwargs = dict()
    kwargs['loadMol'] = _parseFlag(params.get('loadMol', True))
    kwargs['useRDKitChemistry'] = _parseFlag(params.get('useRDKitChemistry', False))
    return _exclude(data, **kwargs)


#-----------------------------------------------------------------------------------------------------------------------

@app.route('/exclude', method=['OPTIONS', 'POST'], name='exclude')
def exclude():
    """
Check if the structure will be excluded from ChEMBL.
Examples and documentation: []()
CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @exclude.mol ${BEAKER_ROOT_URL}exclude
    curl -X POST -F "file=@exclude.mol" ${BEAKER_ROOT_URL}exclude

    """
    data = list(request.files.values())[0].file.read() if len(request.files) else request.body.read()
    return excludeView(data, request.params)

#-----------------------------------------------------------------------------------------------------------------------