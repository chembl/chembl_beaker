from beaker import app
from bottle import request, response
from beaker.core_apps.calculations.impl import  _addHs, _removeHs
from beaker.utils.io import _parseFlag

# ----------------------------------------------------------------------------------------------------------------------


def addHsView(data, params):

    kwargs = dict()
    kwargs['explicitOnly'] = _parseFlag(params.get('explicitOnly', False))
    kwargs['addCoords'] = _parseFlag(params.get('addCoords', False))

    return _addHs(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/addHs', method=['OPTIONS', 'POST'], name="addHs")
def addHs():
    """
Adds hydrogens to the graph of a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @addHs.mol ${BEAKER_ROOT_URL}addHs
    curl -X POST -F "file=@addHs.mol" -F "addCoords=1" ${BEAKER_ROOT_URL}addHs
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return addHsView(data, request.params)

# ----------------------------------------------------------------------------------------------------------------------


def removeHsView(data, params):

    kwargs = dict()
    kwargs['implicitOnly'] = _parseFlag(params.get('implicitOnly', False))
    return _removeHs(data, **kwargs)

# ----------------------------------------------------------------------------------------------------------------------


@app.route('/removeHs', method=['OPTIONS', 'POST'], name="removeHs")
def removeHs():
    """
Removes any hydrogens from the graph of a molecule. CTAB is either single molfile or SDF file.
cURL examples:

    curl -X POST --data-binary @removeHs.mol ${BEAKER_ROOT_URL}removeHs
    curl -X POST -F "file=@removeHs.mol" -F "implicitOnly=1" ${BEAKER_ROOT_URL}removeHs
    """

    data = request.files.values()[0].file.read() if len(request.files) else request.body.read()
    return removeHsView(data, request.params)
